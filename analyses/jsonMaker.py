import supy,steps,samples,calculables

class jsonMaker(supy.analysis) :
    def parameters(self) :
        jwPrompt = calculables.other.jsonWeight("cert/Cert_160404-178677_7TeV_PromptReco_Collisions11_JSON.sub.txt")
        jwMay = calculables.other.jsonWeight("cert/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt")
        jwAug = calculables.other.jsonWeight("cert/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt")

        jw2012 = calculables.other.jsonWeight("cert/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt")
        jw2012_v2 = calculables.other.jsonWeight("cert/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt")
        
        group = self.vary()

        group['SingleMu'] = [(["SingleMu.Run2012A-13Jul2012-v1.AOD.job358",
                               "SingleMu.Run2012A-13Jul2012-v1.AOD.job467",
                               
                               "SingleMu.Run2012A-recover-06Aug2012-v1.AOD.job359",
                               "SingleMu.Run2012A-recover-06Aug2012-v1.AOD.job477", 
                               
                               "SingleMu.Run2012B-13Jul2012-v1.AOD.job358",
                               "SingleMu.Run2012B-13Jul2012-v1.AOD.job461",                                                                
                               "SingleMu.Run2012C-24Aug2012-v1.AOD.job361",
                               "SingleMu.Run2012C-24Aug2012-v1.AOD.job470",
                               
                               "SingleMu.Run2012C-PromptReco-v2.AOD.job360",
                               "SingleMu.Run2012C-PromptReco-v2.AOD.job747",
                               "SingleMu.Run2012D-PromptReco-v1.AOD.job508",
                               "SingleMu.Run2012D-PromptReco-v1.AOD.job525"
                                ], jw2012_v2)]

#        group['SingleEl'] = [(['SingleElectron.Run2012A-PromptReco-v1.AOD.job229',
#                               'SingleElectron.Run2012B-PromptReco-v1.AOD.job228',
#                               'SingleElectron.Run2012B-PromptReco-v1.AOD.job238',
#                               ], jw2012)] 
#
#        group['Photon'] = [(["Photon.Run2012A-PromptReco-v1.AOD.job228",
#                             "SinglePhoton.Run2012B-PromptReco-v1.AOD.job229",
#                             "SinglePhoton.Run2012B-PromptReco-v1.AOD.job238",
#                             ], jw2012)]
#
#        group['Mumu'] = [(['DoubleMu.Run2012A-PromptReco-v1.AOD.job229',
#                           'DoubleMu.Run2012B-PromptReco-v1.AOD.job228',
#                           'DoubleMu.Run2012B-PromptReco-v1.AOD.job239',
#                           ], jw2012)]
#
#        group['HT1'] = [(['HT.Run2012A-PromptReco-v1.AOD.job229',
#                          'HTMHT.Run2012B-PromptReco-v1.AOD.job228',
#                          'HTMHT.Run2012B-PromptReco-v1.AOD.job238',
#                          'JetHT.Run2012B-PromptReco-v1.AOD.job228',
#                          'JetHT.Run2012B-PromptReco-v1.AOD.job238',
#                          ], jw2012)]
#        
#        group['HT2'] = [('HT.Run2011A-May10ReReco-v1.AOD.Darren1',jwMay),
#                        ('HT.Run2011A-05Aug2011-v1.AOD.Bryn1',jwAug),
#                        (['HT.Run2011A-PromptReco-v4.AOD.Bryn1',
#                          'HT.Run2011A-PromptReco-v6.AOD.Bryn1',
#                          'HT.Run2011B-PromptReco-v1.AOD.Bryn1',
#                          'HT.Run2011B-PromptReco-v1.AOD.Bryn2',
#                          'HT.Run2011B-PromptReco-v1.AOD.Bryn3',
#                          ], jwPrompt)]
#
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
