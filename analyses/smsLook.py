import supy,steps,calculables,samples,ROOT as r

class smsLook(supy.analysis) :
    def parameters(self) :
        return {}

    def listOfCalculables(self, params) :
        return supy.calculables.zeroArgs(calculables)+supy.calculables.zeroArgs(supy.calculables)
    
    def listOfSteps(self, params) :
        return [
            supy.steps.printer.progressPrinter(),
            supy.steps.other.collector(["SimpModelScanmGL","SimpModelScanmLSP"]),
            #steps.other.smsMedianHistogrammer(_jet),
            ]#+[supy.steps.histos.generic(("SimpModelScanmGL","SimpModelScanmLSP"),(101,38),(i,i),(2020+i,760+i),title="%d;m_{0};m_{1/2};events/bin"%i) for i in range(20)]
    
    def listOfSampleDictionaries(self) :
        sampleDict = supy.samples.SampleHolder()
        sampleDict.add("t1_1000_50", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_50/t1_1000_50.root"]', lumi = 1.1e3)
        sampleDict.add("t1_1000_600", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_600/t1_1000_600.root"]', lumi = 1.1e3)
        sampleDict.add("t1_400_300", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim400_300/t1_400_300.root"]', lumi = 1.1e3)        
        sampleDict.add("t1_3_points", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim/sms_3_points.root"]', lumi = 1.1e3)

        return [sampleDict,samples.susy17]

    def listOfSamples(self,params) :
        from supy.samples import specify
        return supy.samples.specify(names = ["T2tt_8.job351"], nFilesMax = 1)
   
    def conclude(self,pars) :
        org = self.organizer(pars)
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          doLog = False,
                          #pegMinimum = 0.1,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          ).plotAll()
