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
                calculables.gen.genIndices( pdgs = [-16,-14,-12,12,14,16], label = "Status3Nu", status = [3]),
                calculables.gen.genIndicesPtSorted(label = "Status3b"),
                calculables.gen.genIndicesPtSorted(label = "Status3t"),
                calculables.gen.genIndicesPtSorted(label = "Status3W"),
                calculables.gen.genSumPt(indexLabels = ["Status3b", "Status3W"]),
                calculables.gen.genSumPt(indexLabels = ["Status3b", "Status3WDaughters"]),
                calculables.gen.genIndicesHardPartons(),

                calculables.gen.JetIndices(("ak5Gen", ""), ptMin = 10.0, etaMax = 3.0),
                ]
        return out

    def listOfSampleDictionaries(self) :
        return [samples.top17]

    def listOfSamples(self, params) :
        return supy.samples.specify(names = "ttz_8_mg.job269_1")

    def triggerFilters(self, thresh = tuple()) :
        "/online/collisions/2012/7e33/v4.0/HLT/V11"
        out = []

        if thresh==(80, 80, 80, 80) :
            "HLT_QuadJet80_v6"
            out = [supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 3, min = 80.0)]

        elif thresh==(80, 80, 60, 60, 20, 20) :
            "HLT_DiJet80_DiJet60_DiJet20_v5"
            out = [supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 1, min = 80.0),
                   supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 3, min = 60.0),
                   supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 5, min = 20.0),
                   ]

        elif thresh==(60, 60, 60, 60, 20, 20) :
            "HLT_QuadJet60_DiJet20_v5"
            out = [supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 3, min = 60.0),
                   supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 5, min = 20.0),
                   ]

        elif thresh==(45, 45, 45, 45, 45, 45) :
            "HLT_SixJet45_v6"
            out = [supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 5, min = 45.0),
                   ]

        elif thresh==tuple([30]*8) :
            "HLT_EightJet30_eta3p0_v5"
            out = [supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 7, min = 30.0),
                   ]

        return out

    def genFilters(self) :
        return [supy.steps.filters.multiplicity("genIndicesStatus3Z", min = 1, max = 1),
                supy.steps.filters.multiplicity("genIndicesStatus3t", min = 2, max = 2),
                supy.steps.filters.multiplicity("genIndicesStatus3W", min = 2, max = 2),
                supy.steps.filters.multiplicity("genIndicesStatus3b", min = 4, max = 4),
                supy.steps.filters.multiplicity("genIndicesStatus3WDaughters", min = 4, max = 4),
                supy.steps.filters.multiplicity("genIndicesStatus3Nu", max = 0),
                ]

    def scaleLook(self) :
        return [supy.steps.histos.value("genRootSHat", 20, 0.0, 8000.0, xtitle = "#sqrt{#hat{s}} (GeV)"),
                supy.steps.histos.value("genSumPt_Status3b_Status3W",          20, 0.0, 1000.0),
                supy.steps.histos.value("genSumPt_Status3b_Status3WDaughters", 20, 0.0, 1000.0),
                ]

    def printers(self) :
        return [supy.steps.printer.progressPrinter(),
                #steps.gen.particlePrinter(minPt = -1.0, minStatus = 3),
                steps.gen.genJetPrinter(cs = ("genak5", "")),
                ]

    def listOfSteps(self, params) :
        outList = []

        outList += self.genFilters()
        #outList += self.printers()
        #outList += self.scaleLook()
        outList += self.triggerFilters(thresh = (80, 80, 60, 60, 20, 20))

        outList += [
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
