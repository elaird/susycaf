import supy,steps,samples,calculables

class jsonMaker16(supy.analysis) :
    def parameters(self) :

        group = self.vary()

        group['SingleMu'] = [(["SingleMu.2011A",
                               "SingleMu.2011B",
                               ], [])]

        group['MuHad'] = [(["MuHad.2011A",
                            "MuHad.2011B",
                            ],[])]

        group['SingleEl'] = [(["SingleEl.2011A",
                               "SingleEl.2011B",
                               ], [])] 

        group['EleHad'] = [(["EleHad.2011A",
                             "EleHad.2011B",
                             ], [])]

        return {'group':group}

    def listOfSteps(self,pars) :
        return [ supy.steps.printer.progressPrinter(2,300),
                 steps.other.jsonMaker(pixelLumi = True),
                 ]

    def listOfCalculables(self,pars) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,pars) :
        return sum([supy.samples.specify(names = samps, weights = jw) for samps,jw in pars['group']],[])

    def listOfSampleDictionaries(self) :
        return [samples.electron16, samples.muon16]

    def mainTree(self) :
        return ("lumiTree","tree")

    def otherTreesToKeepWhenSkimming(self) :
        return []
