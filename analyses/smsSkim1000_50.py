import supy,calculables,steps,samples


class smsSkim1000_50(supy.analysis) :

    def listOfSteps(self, params) :
        stepList = [
            supy.steps.printer.progressPrinter(2,300),
            supy.steps.filters.value("SimpModelScanmGL", min = 1000, max = 1000),
            supy.steps.filters.value("SimpModelScanmLSP", min = 50, max = 50),
            #supy.steps.filters.multiplicity("genIndicesStatus3Z", min = 1, max = 1),
            #supy.steps.histos.mass("genP4", 100, 0, 200, index = 0, indices = "genIndicesStatus3Z", xtitle="Z_{m} (GeV)"),
            #supy.steps.filters.mass("genP4", index = 0, indices = "genIndicesStatus3Z", min = 65, max = 110),
            #supy.steps.histos.value("xcak5JetSumEtPat",100,0,500,xtitle="H_{T} (GeV)"),
            #steps.jet.htSelector(nameList(params["recoAlgos"], "jet"), 250.0),
            supy.steps.other.skimmer(),
            ]

        return stepList

    def listOfCalculables(self, params) :
        outList = supy.calculables.zeroArgs(supy.calculables)
        return outList
    
    def listOfSamples(self, params) :
        from supy.samples import specify
        
        return specify( #nFilesMax = 1, nEventsMax = 2000,
                       names = ["t1.yos"])

    def listOfSampleDictionaries(self) :
        return [samples.mc]

    
    def conclude(self, config) :
        org = self.organizer(config)
        supy.utils.io.printSkimResults(org)

        supy.plotter(org,
                     pdfFileName = self.pdfFileName(org.tag),
                     blackList = ["lumiHisto","xsHisto","nJobsHisto",],
                     ).plotAll()

