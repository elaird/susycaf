import supy,steps,samples,calculables

class jsonMaker(supy.analysis) :
    def parameters(self) :
        jw2012 = calculables.other.jsonWeight("cert/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt")       
        group = self.vary()

        group['SingleMu'] = [(["SingleMu.Run2012%s-22Jan2013"
                               % era for era in ["A","B","C","D"]],
                              jw2012)]

        group['Photon']=[(["Photon.Run2012A-22Jan2013",
                           "SinglePhotonParked.Run2012B-22Jan2013",
                           "SinglePhotonParked.Run2012C-22Jan2013",
                           "SinglePhotonParked.Run2012D-22Jan2013"],
                          jw2012)]

        group['HT'] = [(["HT.Run2012A-22Jan2013"] +
                        ["HTMHTParked.Run2012%s-22Jan2013"
                         % era for era in ["B","C","D"]],
                        jw2012)]

        return {'group':group}
    
    def listOfSteps(self,pars) :
        return [ supy.steps.printer.progressPrinter(2,300),
                 steps.other.jsonMaker(pixelLumi = False),
                 ]

    def listOfCalculables(self,pars) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,pars) :
        return sum([supy.samples.specify(names = samps, weights = jw) for samps,jw in pars['group']],[])

    def listOfSampleDictionaries(self) :
        return [samples.ht17, samples.jetmet17, samples.muon17, samples.photon17, samples.electron17, samples.mumu17]

    def mainTree(self) :
        return ("lumiTree","tree")

    def otherTreesToKeepWhenSkimming(self) :
        return []
