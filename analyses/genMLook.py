import supy,steps,calculables,samples
import ROOT as r

class genMLook(supy.analysis) :

    def listOfCalculables(self, params) :
        out = supy.calculables.zeroArgs(supy.calculables) + supy.calculables.zeroArgs(calculables)
        out += [calculables.gen.genIndices( pdgs = [23], label = "Status3Z", status = [3]),
                ]
        return out

    def listOfSampleDictionaries(self) :
        return [samples.top17]

    def listOfSamples(self, params) :
        return supy.samples.specify(names = "ttz_8_mg.job269_1")

    def listOfSteps(self, params) :
        outList=[
            supy.steps.printer.progressPrinter(),
            steps.gen.particlePrinter(minPt = -1.0, minStatus = 3),
            supy.steps.filters.multiplicity("genIndicesStatus3Z", min = 1, max = 1),
            supy.steps.histos.mass("genP4", 100, 0.0, 200.0, indices = "genIndicesStatus3Z", index = 0, xtitle = "gen Z0"),
            #steps.gen.genMassHistogrammer(),
            #steps.gen.genSHatHistogrammer(),
            ]
        return outList

    def conclude(self,pars) :
        org = self.organizer(pars)
        org.scale(1000.0)
        
        supy.plotter(org,
                     pdfFileName = self.pdfFileName(""),
                     blackList = ["lumiHisto","xsHisto","nJobsHisto",],
                     rowColors = [r.kBlack, r.kViolet+4],
                     pegMinimum = 0.1,
                     ).plotAll()
