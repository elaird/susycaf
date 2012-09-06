import supy,steps,calculables,samples
import ROOT as r

class genTopLook(supy.analysis) :

    def listOfCalculables(self, params) :
        out = supy.calculables.zeroArgs(supy.calculables) + supy.calculables.zeroArgs(calculables)
        out += [calculables.gen.genIndices( pdgs = [-5,5], label = "b", status = [3]),
                calculables.gen.genIndices( pdgs = [-24,24], label = "W", status = [3]),
                calculables.gen.genIndices( pdgs = [-6,6], label = "t", status = [3]),
                calculables.gen.genIndices( pdgs = range(-16,-10)+range(11,17), label = "lepton", status = [3]),

                calculables.gen.genIndices( pdgs = [ 5], label = "b+_t+Daughters", status = [3], motherPdgs = [  6]),
                calculables.gen.genIndices( pdgs = [  ], label = "W+Daughters",    status = [3], motherPdgs = [ 24]),
                calculables.gen.genIndices( pdgs = [-5], label = "b-_t-Daughters", status = [3], motherPdgs = [ -6]),
                calculables.gen.genIndices( pdgs = [  ], label = "W-Daughters",    status = [3], motherPdgs = [-24]),
                ]
        return out

    def listOfSampleDictionaries(self) :
        return [samples.top17]

    def listOfSamples(self, params) :
        return supy.samples.specify(names = "tt_1", nFilesMax = 1)

    def listOfSteps(self, params) :
        outList=[
            supy.steps.printer.progressPrinter(),
            supy.steps.filters.multiplicity("genIndicest", min = 2, max = 2),
            supy.steps.filters.multiplicity("genIndicesW", min = 2, max = 2),
            supy.steps.filters.multiplicity("genIndicesb", min = 2, max = 2),
            supy.steps.filters.multiplicity("genIndiceslepton", max = 0),
            supy.steps.filters.multiplicity("genIndicesb+_t+Daughters", min = 1, max = 1),
            supy.steps.filters.multiplicity("genIndicesb-_t-Daughters", min = 1, max = 1),
            supy.steps.filters.multiplicity("genIndicesW+Daughters", min = 2, max = 2),
            supy.steps.filters.multiplicity("genIndicesW-Daughters", min = 2, max = 2),

            #steps.gen.particlePrinter(minPt = -1.0, minStatus = 3),
            steps.gen.bqDotHistogrammer(),
            ]
        return outList

    def conclude(self,pars) :
        org = self.organizer(pars)
        org.scale(1.0e3)

        supy.plotter(org,
                     pdfFileName = self.pdfFileName(""),
                     blackList = ["lumiHisto","xsHisto","nJobsHisto",],
                     rowColors = [r.kBlack, r.kViolet+4],
                     pegMinimum = 0.1,
                     doLog = False,
                     ).plotAll()
