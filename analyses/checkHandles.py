import supy
import steps,samples

class checkHandles(supy.analysis) :

    def listOfSteps(self,params) :
        stepList=[steps.printer.progressPrinter(2,300),
                  #steps.Other.handleChecker(),
                  supy.steps.histos.iterHistogrammer("ak5JetCorrectedP4Pat", 100, 0.0, 100.0, title=";ak5 calo jet corrected p_{T} (GeV);jets / bin", funcString="lambda x:x.pt()"),
                  supy.steps.histos.iterHistogrammer("ak5JetPFCorrectedP4Pat", 100, 0.0, 100.0, title=";ak5 PF jet corrected p_{T} (GeV);jets / bin", funcString="lambda x:x.pt()"),
                  supy.steps.histos.iterHistogrammer("electronP4Pat", 100, 0.0, 100.0, title=";electron p_{T} (GeV);electrons / bin", funcString="lambda x:x.pt()"),
                  supy.steps.histos.iterHistogrammer("photonP4Pat", 100, 0.0, 100.0, title=";photon p_{T} (GeV);photons / bin", funcString="lambda x:x.pt()"),
                  supy.steps.histos.iterHistogrammer("muonP4Pat", 100, 0.0, 100.0, title=";muon p_{T} (GeV);muons / bin", funcString="lambda x:x.pt()"),
                  supy.steps.histos.iterHistogrammer("tauP4Pat", 100, 0.0, 100.0, title=";tau p_{T} (GeV);taus / bin", funcString="lambda x:x.pt()"),
                  supy.steps.histos.iterHistogrammer("genak5GenJetsP4", 100, 0.0, 100.0, title=";gen. ak5 jet corrected p_{T} (GeV);jets / bin", funcString="lambda x:x.pt()"),
                  supy.steps.histos.iterHistogrammer("genStatus1P4", 100, 0.0, 100.0, title=";status 1 gen. particle p_{T} (GeV);jets / bin", funcString="lambda x:x.pt()"),
                  ]
        return stepList

    def listOfCalculables(self,params) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,params) :
        return ( supy.samples.specify(names = "41testData") + 
                 supy.samples.specify(names = "41testMc", color = 2) )

    def listOfSampleDictionaries(self) :
        return [samples.jetmet]

    def conclude(self,pars) :
        org = self.organizer(pars)
        supy.plotter(org, psFileName = self.psFileName( org.tag ), doLog = False,
                     blackList = ["lumiHisto","xsHisto","nJobsHisto"]).plotAll()
