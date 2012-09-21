import supy,steps,calculables,samples,ROOT as r

class smsLook(supy.analysis) :
    def parameters(self) :
        return {}

    def listOfCalculables(self, params) :
        out  = []
        out += supy.calculables.zeroArgs(calculables)
        out += supy.calculables.zeroArgs(supy.calculables)
        out += [calculables.gen.genIndices( pdgs = [-16,-14,-12,12,14,16], label = "Status3Nu", status = [3]),
                calculables.gen.genIndices( pdgs = [-6,6], label = "Status3t", status = [3]),
                calculables.gen.SumP4(indices = "genIndicesStatus3t"),
                ]
        return out
                
    def stepsPrepare(self, params) :
        return [supy.steps.other.collector(["susyScanmGL","susyScanmLSP"]),
                supy.steps.filters.value('susyScanmGL', min = 499, max = 501),
                supy.steps.filters.value('susyScanmLSP', min = 99, max = 101),
                supy.steps.other.skimmmer(),
                ]

    def listOfSteps(self, params) :
        return [
            supy.steps.printer.progressPrinter(),
            #self.stepsPrepare(params),
            
            supy.steps.histos.pt('genmetP4True', 100, 0.0, 500.0, xtitle = "gen. MET (GeV)"),
            supy.steps.histos.multiplicity('genIndicesStatus3Nu'),
            supy.steps.filters.multiplicity('genIndicesStatus3Nu', max = 0),
            supy.steps.histos.pt('genmetP4True', 100, 0.0, 500.0, xtitle = "gen. MET (GeV)"),

            supy.steps.histos.multiplicity('genIndicesStatus3t'),
            supy.steps.histos.pt('SumP4genIndicesStatus3t', 100, 0.0, 500.0, xtitle = "(t+t system)"),

            supy.steps.histos.histogrammer(("genmetP4True", "SumP4genIndicesStatus3t"), (100, 100), (0.0, 0.0), (1000.0, 1000.0),
                                           title = ";gen. MET (GeV);(t+t system) p_{T} (GeV);events / bin",
                                           funcString = "lambda x:(x[0].pt(),x[1].pt())"),

            #steps.other.smsMedianHistogrammer(_jet),
            ]#+[supy.steps.histos.generic(("SimpModelScanmGL","SimpModelScanmLSP"),(101,38),(i,i),(2020+i,760+i),title="%d;m_{0};m_{1/2};events/bin"%i) for i in range(20)]
    
    def listOfSampleDictionaries(self) :
        sampleDict = supy.samples.SampleHolder()
        sampleDict.add("t1_1000_50", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_50/t1_1000_50.root"]', lumi = 1.1e3)
        sampleDict.add("t1_1000_600", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_600/t1_1000_600.root"]', lumi = 1.1e3)
        sampleDict.add("t1_400_300", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim400_300/t1_400_300.root"]', lumi = 1.1e3)        
        sampleDict.add("t1_3_points", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim/sms_3_points.root"]', lumi = 1.1e3)

        return [sampleDict,samples.susy17,samples.top17]

    def listOfSamples(self,params) :
        from supy.samples import specify
        return supy.samples.specify(names = ["T2tt_500_100"])+supy.samples.specify(names = ["tt_8_mg.job315_1"])

    def conclude(self,pars) :
        org = self.organizer(pars)

        def md(x, y) :
            x.update(y)
            return x
        mcOps = {"markerStyle":1, "lineWidth":1, "goptions":"ehist"}
        org.mergeSamples(targetSpec = md({"name":"stop-stop (500,100)", "color":r.kRed}, mcOps), allWithPrefix = "T2tt_500_100")
        org.mergeSamples(targetSpec = md({"name":"t-t", "color":r.kBlue}, mcOps), allWithPrefix = "tt")

        org.scale(20.0e3)
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          #doLog = False,
                          pegMinimum = 0.1,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          ).plotAll()
