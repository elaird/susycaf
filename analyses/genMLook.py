import supy,steps,calculables,samples
import ROOT as r

class genMLook(supy.analysis) :

    def listOfCalculables(self, params) :
        out = supy.calculables.zeroArgs(supy.calculables) + supy.calculables.zeroArgs(calculables)
        out += [calculables.gen.genIndices( pdgs = [23], label = "Status3Z", status = [3]),
                calculables.gen.genIndices( pdgs = [-5,5], label = "Status3b", status = [3]),
                calculables.gen.genIndices( pdgs = [-24,24], label = "Status3W", status = [3]),
                calculables.gen.genIndices( pdgs = [], label = "Status3WDaughters", status = [3], motherPdgs = [-24, 24]),
                calculables.gen.genIndices( pdgs = [-6,6], label = "Status3t", status = [3]),
                calculables.gen.genIndicesPtSorted(label = "Status3b"),
                calculables.gen.genIndicesPtSorted(label = "Status3t"),
                calculables.gen.genIndicesPtSorted(label = "Status3W"),
                calculables.gen.genSumPt(indexLabels = ["Status3b", "Status3W"]),
                calculables.gen.genSumPt(indexLabels = ["Status3b", "Status3WDaughters"]),
                ]
        return out

    def listOfSampleDictionaries(self) :
        return [samples.top17]

    def listOfSamples(self, params) :
        return supy.samples.specify(names = "ttz_8_mg.job269_1")

    def listOfSteps(self, params) :
        outList=[
            supy.steps.printer.progressPrinter(),
            #steps.gen.particlePrinter(minPt = -1.0, minStatus = 3),
            supy.steps.filters.multiplicity("genIndicesStatus3Z", min = 1, max = 1),
            supy.steps.filters.multiplicity("genIndicesStatus3t", min = 2, max = 2),
            supy.steps.filters.multiplicity("genIndicesStatus3W", min = 2, max = 2),
            supy.steps.filters.multiplicity("genIndicesStatus3b", min = 4, max = 4),
            supy.steps.filters.multiplicity("genIndicesStatus3WDaughters", min = 4, max = 4),

            supy.steps.histos.value("genRootSHat", 20, 0.0, 8000.0, xtitle = "#sqrt{#hat{s}} (GeV)"),
            supy.steps.histos.value("genSumPt_Status3b_Status3W",          20, 0.0, 1000.0),
            supy.steps.histos.value("genSumPt_Status3b_Status3WDaughters", 20, 0.0, 1000.0),

            supy.steps.histos.mass("genP4", 100, 0.0, 200.0, indices = "genIndicesStatus3Z", index = 0, xtitle = "gen Z_{0}"),
            supy.steps.histos.pt  ("genP4",  20, 0.0, 300.0, indices = "genIndicesStatus3Z", index = 0, xtitle = "gen Z_{0}"),

            supy.steps.histos.pt("genP4", 20, 0.0, 300.0, indices = "genIndicesStatus3tPtSorted", index = 0, xtitle = "    leading t quark"),
            supy.steps.histos.pt("genP4", 20, 0.0, 300.0, indices = "genIndicesStatus3tPtSorted", index = 1, xtitle = "sub-leading t quark"),

            supy.steps.histos.pt("genP4", 20, 0.0, 300.0, indices = "genIndicesStatus3WPtSorted", index = 0, xtitle = "    leading W"),
            supy.steps.histos.pt("genP4", 20, 0.0, 300.0, indices = "genIndicesStatus3WPtSorted", index = 1, xtitle = "sub-leading W"),

            supy.steps.histos.pt("genP4", 20, 0.0, 200.0, indices = "genIndicesStatus3bPtSorted", index = 0, xtitle = "    leading b quark"),
            supy.steps.histos.pt("genP4", 20, 0.0, 200.0, indices = "genIndicesStatus3bPtSorted", index = 1, xtitle = "sub-leading b quark"),
            supy.steps.histos.pt("genP4", 20, 0.0, 200.0, indices = "genIndicesStatus3bPtSorted", index = 2, xtitle = "3rd leading b quark"),
            supy.steps.histos.pt("genP4", 20, 0.0, 200.0, indices = "genIndicesStatus3bPtSorted", index = 3, xtitle = "4th leading b quark"),
            #steps.gen.genSHatHistogrammer(),
            ]
        return outList

    def conclude(self,pars) :
        org = self.organizer(pars)
        org.mergeSamples(sources = ["ttz_8_mg.job269_1"], targetSpec = {"name":"TTZ", "markerStyle":20})
        org.scale(10.0e3)

        supy.plotter(org,
                     pdfFileName = self.pdfFileName(""),
                     blackList = ["lumiHisto","xsHisto","nJobsHisto",],
                     rowColors = [r.kBlack, r.kViolet+4],
                     pegMinimum = 0.1,
                     doLog = False,
                     ).plotAll()
