import supy,steps,calculables,samples,ROOT as r

class smsLook(supy.analysis) :
    def parameters(self) :
        return {}

    def listOfCalculables(self, params) :
        out  = []
        out += supy.calculables.zeroArgs(calculables)
        out += supy.calculables.zeroArgs(supy.calculables)
        out += [calculables.gen.Indices(pdgs = [-16,-14,-12,12,14,16], label = "Status3Nu", status = [3]),
                calculables.gen.Indices(pdgs = [-6,6], label = "Status3t", status = [3]),
                calculables.gen.Indices(pdgs = [-24,24], label = "Status3W", status = [3]),
                calculables.gen.Indices(pdgs = [ 6], label = "Status3t+", status = [3]),
                calculables.gen.Indices(pdgs = [-6], label = "Status3t-", status = [3]),
                calculables.gen.Indices(pdgs = [ 24], label = "Status3W+", status = [3]),
                calculables.gen.Indices(pdgs = [-24], label = "Status3W-", status = [3]),
                calculables.gen.Indices(pdgs = [-5,5], label = "Status3b", status = [3]),
                calculables.gen.Indices(pdgs = [-11,11], label = "Status3e", status = [3]),
                calculables.gen.Indices(pdgs = [-13,13], label = "Status3mu", status = [3]),
                calculables.gen.Indices(pdgs = [-15,15], label = "Status3tau", status = [3]),
                calculables.gen.Indices(pdgs = [-15,15], label = "2TauFrom3Tau", status = [2], motherIndices = "IndicesStatus3tau"),

                calculables.gen.Indices(pdgs = [-11,11], label = "Status1eFromTau",
                                        status = [1], motherIndices =  "Indices2TauFrom3Tau"),
                calculables.gen.Indices(pdgs = [-13,13], label = "Status1muFromTau",
                                        status = [1], motherIndices = "Indices2TauFrom3Tau"),

                calculables.gen.Indices(pdgs = [], label = "Status3t+Daughters", status = [3], motherPdgs = [ 6]),
                calculables.gen.Indices(pdgs = [], label = "Status3t-Daughters", status = [3], motherPdgs = [-6]),
                calculables.gen.Indices(pdgs = [], label = "Status3W+Daughters", status = [3], motherPdgs = [ 24]),
                calculables.gen.Indices(pdgs = [], label = "Status3W-Daughters", status = [3], motherPdgs = [-24]),

                calculables.gen.IndicesPtSorted(label = "Status3t"),
                calculables.gen.IndicesPtSorted(label = "Status3W"),
                calculables.gen.IndicesPtSorted(label = "Status3b"),

                calculables.gen.KineIndices(indices = "IndicesStatus3b", ptMin = 30.0, etaMax = 2.2),
                calculables.gen.KineIndices(indices = "IndicesStatus3e", ptMin = 10.0, etaMax = 2.3),
                calculables.gen.KineIndices(indices = "IndicesStatus3mu", ptMin = 10.0, etaMax = 2.3),
                calculables.gen.KineIndices(indices = "IndicesStatus1eFromTau", ptMin = 10.0, etaMax = 2.3),
                calculables.gen.KineIndices(indices = "IndicesStatus1muFromTau", ptMin = 10.0, etaMax = 2.3),

                supy.calculables.other.pt("genP4", indices = "IndicesStatus3t+", index = 0),
                supy.calculables.other.pt("genP4", indices = "IndicesStatus3t-", index = 0),

                calculables.gen.SumP4(indices = "IndicesStatus3t"),
                calculables.gen.DeltaR(indices = "IndicesStatus3t+Daughters"),
                calculables.gen.DeltaR(indices = "IndicesStatus3t-Daughters"),
                calculables.gen.DeltaR(indices = "IndicesStatus3W+Daughters"),
                calculables.gen.DeltaR(indices = "IndicesStatus3W-Daughters"),

                calculables.gen.MinDeltaPhiMet(indices = "IndicesStatus3b", met = "genmetP4True"),
                calculables.gen.MinDeltaPhiMet(indices = "IndicesStatus3W+Daughters", met = "genmetP4True"),
                calculables.gen.MinDeltaPhiMet(indices = "IndicesStatus3W-Daughters", met = "genmetP4True"),

                calculables.gen.JetIndices(("ak5Gen", ""), ptMin = 20.0, etaMax = 3.0),
                ]
        return out
                
    def stepsPrepare(self, params) :
        return [supy.steps.other.collector(["susyScanmGL","susyScanmLSP"]),
                supy.steps.filters.value('susyScanmGL', min = 499, max = 501),
                supy.steps.filters.value('susyScanmLSP', min = 299, max = 301),
                supy.steps.other.skimmer(),
                ]

    def triggerFilters(self, thresh = tuple(), offline = False) :
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
            if offline :
                out += [
                    supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 3, min = 80.0),
                    supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 5, min = 30.0),
                    ]

        elif thresh==(45, 45, 45, 45, 45, 45) :
            "HLT_SixJet45_v6"
            out = [supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 5, min = 45.0),
                   ]

        elif thresh==tuple([30]*8) :
            "HLT_EightJet30_eta3p0_v5"
            out = [supy.steps.filters.pt("genak5GenJetsP4", indices = "ak5GenJetIndices", index = 7, min = 30.0),
                   ]

        elif thresh==(150,) :
            "HLT_PFMET150_v6"
            out = [supy.steps.filters.pt("genmetP4True", min = 150.0),]
            if offline :
                out += [supy.steps.filters.pt("genmetP4True", min = 180.0),]

        if out : out.insert(0, supy.steps.filters.label('trigger requirements'))
        return out

    def genFilters(self) :
        return [
            #supy.steps.histos.pt('genmetP4True', 100, 0.0, 500.0, xtitle = "gen. MET (GeV)"),
            #supy.steps.histos.multiplicity('IndicesStatus3Nu'),
            #supy.steps.filters.multiplicity('IndicesStatus3Nu', max = 0),
            supy.steps.filters.pt('genmetP4True', min = 100.0),
            #supy.steps.histos.multiplicity('IndicesStatus3t'),
            #supy.steps.other.skimmer(),
            ]
    def met(self) :
        return [supy.steps.histos.pt('genmetP4True', 100, 0.0, 500.0, xtitle = "gen. MET (GeV)"),]

    def b(self) :
        return [supy.steps.filters.label('b requirements'),
                supy.steps.filters.multiplicity("KineIndicesStatus3b", min = 2, max = 2),
                ]

    def jet(self) :
        return [supy.steps.filters.label('jet requirements'),
                #supy.steps.histos.multiplicity('ak5GenJetIndices', max = 15),
                supy.steps.filters.multiplicity('ak5GenJetIndices', min = 6),
                supy.steps.filters.multiplicity('ak5GenJetIndices', max = 8),
                ]

    def lepton(self) :
        return [
            supy.steps.filters.label('lepton vetoes'),
            #supy.steps.histos.multiplicity("KineIndicesStatus3e"),
            #supy.steps.histos.multiplicity("KineIndicesStatus3mu"),
            #supy.steps.histos.pt("genP4", 100, 0.0, 100.0, indices = "KineIndicesStatus3e", index = 0),
            #supy.steps.histos.pt("genP4", 100, 0.0, 100.0, indices = "KineIndicesStatus3mu", index = 0),
            supy.steps.filters.multiplicity("KineIndicesStatus3e", max = 0),
            supy.steps.filters.multiplicity("KineIndicesStatus3mu", max = 0),
            supy.steps.filters.multiplicity("KineIndicesStatus1eFromTau", max = 0),
            supy.steps.filters.multiplicity("KineIndicesStatus1muFromTau", max = 0),
            ]

    def ptPlots(self) :
        return [
            supy.steps.histos.pt("genP4", 40, 0.0, 400.0, indices = "IndicesStatus3tPtSorted", index = 0, xtitle = "    leading t quark"),
            supy.steps.histos.pt("genP4", 40, 0.0, 400.0, indices = "IndicesStatus3tPtSorted", index = 1, xtitle = "sub-leading t quark"),

            supy.steps.histos.pt("genP4", 40, 0.0, 400.0, indices = "IndicesStatus3WPtSorted", index = 0, xtitle = "    leading W"),
            supy.steps.histos.pt("genP4", 40, 0.0, 400.0, indices = "IndicesStatus3WPtSorted", index = 1, xtitle = "sub-leading W"),

            supy.steps.histos.pt("genP4", 20, 0.0, 200.0, indices = "IndicesStatus3bPtSorted", index = 0, xtitle = "    leading b quark"),
            supy.steps.histos.pt("genP4", 20, 0.0, 200.0, indices = "IndicesStatus3bPtSorted", index = 1, xtitle = "sub-leading b quark"),

            supy.steps.histos.absEta("genP4", 20, 0.0, 3.0, indices = "IndicesStatus3bPtSorted", index = 0, xtitle = "    leading b quark"),
            supy.steps.histos.absEta("genP4", 20, 0.0, 3.0, indices = "IndicesStatus3bPtSorted", index = 1, xtitle = "sub-leading b quark"),

            supy.steps.histos.pt("genak5GenJetsP4", 20, 0.0, 200.0, indices = "ak5GenJetIndices", index = 3, xtitle = "4th gen. jet"),
            supy.steps.histos.pt("genak5GenJetsP4", 20, 0.0, 200.0, indices = "ak5GenJetIndices", index = 5, xtitle = "6th gen. jet"),
            ]

    def deltaRPlots(self) :
        return [
            supy.steps.histos.value('IndicesStatus3t+DaughtersDeltaR', 100, 0.0, 5.0, xtitle = "#DeltaR(b,W) from t"),
            supy.steps.histos.value('IndicesStatus3t-DaughtersDeltaR', 100, 0.0, 5.0, xtitle = "#DeltaR(b,W) from #bar{t}"),
            supy.steps.histos.value('IndicesStatus3W+DaughtersDeltaR', 100, 0.0, 5.0, xtitle = "#DeltaR(q,q') from W+"),
            supy.steps.histos.value('IndicesStatus3W-DaughtersDeltaR', 100, 0.0, 5.0, xtitle = "#DeltaR(q,q') from W-"),

            supy.steps.histos.histogrammer(("genP4.pt0IndicesStatus3t-", "IndicesStatus3t-DaughtersDeltaR"),
                                           (100, 100), (0.0, 0.0), (500.0, 5.0),
                                           title = ";p_{T} of #bar{t} (GeV);#DeltaR(b,W) from #bar{t};events / bin"),
            supy.steps.histos.histogrammer(("genP4.pt0IndicesStatus3t+", "IndicesStatus3t+DaughtersDeltaR"),
                                           (100, 100), (0.0, 0.0), (500.0, 5.0),
                                           title = ";p_{T} of t (GeV);#DeltaR(b,W) from t;events / bin"),
            ]

    def deltaPhi(self) :
        return [
            supy.steps.filters.label('Dphi requirements'),
            #supy.steps.histos.histogrammer("IndicesStatus3bMinDeltaPhiMetgenmetP4True", 100, 0.0, r.TMath.Pi(),
            #                               title = ";min. |#Delta#phi(b,MET)|;events / bin"),
            #supy.steps.histos.histogrammer(("genmetP4True", "IndicesStatus3bMinDeltaPhiMetgenmetP4True"),
            #                               (100, 100), (0.0, 0.0), (300.0, r.TMath.Pi()),
            #                               title = ";MET;min. |#Delta#phi(b,MET)|;events / bin", funcString = "lambda x:(x[0].pt(),x[1])"),
            supy.steps.filters.value("IndicesStatus3bMinDeltaPhiMetgenmetP4True", min = 0.2),
            supy.steps.filters.value("IndicesStatus3W-DaughtersMinDeltaPhiMetgenmetP4True", min = 0.2),
            supy.steps.filters.value("IndicesStatus3W+DaughtersMinDeltaPhiMetgenmetP4True", min = 0.2),
            ]

    def ttPlots(self) :
        return [
            supy.steps.histos.pt('SumP4IndicesStatus3t', 100, 0.0, 500.0, xtitle = "(t+t system)"),
            supy.steps.histos.histogrammer(("genmetP4True", "SumP4IndicesStatus3t"), (100, 100), (0.0, 0.0), (1000.0, 1000.0),
                                           title = ";gen. MET (GeV);(t+t system) p_{T} (GeV);events / bin",
                                           funcString = "lambda x:(x[0].pt(),x[1].pt())"),
            ]

    def progress(self) :
        return [supy.steps.printer.progressPrinter()]

    def printer(self) :
        return [
            steps.printer.eventPrinter(),
            steps.printer.metPrinter(collections = ["genmetP4True"]),
            steps.gen.particlePrinter(),
            ]

    def listOfSteps(self, params) :
        return (
            #self.stepsPrepare(params) +
            self.progress() +
            self.genFilters() +
            self.met() +
            self.triggerFilters(thresh = (150,), offline = True) +
            #self.triggerFilters(thresh = (60, 60, 60, 60, 20, 20), offline = True) +
            self.jet() +
            self.b() +
            self.lepton() +
            #self.deltaPhi() +
            [])+[
            supy.steps.filters.label('misc plots'),
            ]+(
            self.printer() +
            #self.ttPlots() +
            self.met() +
            #self.ptPlots() +
            #self.deltaRPlots() +
            [])

    def listOfSampleDictionaries(self) :
        sampleDict = supy.samples.SampleHolder()
        sampleDict.add("t1_1000_50", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_50/t1_1000_50.root"]', lumi = 1.1e3)
        sampleDict.add("t1_1000_600", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_600/t1_1000_600.root"]', lumi = 1.1e3)
        sampleDict.add("t1_400_300", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim400_300/t1_400_300.root"]', lumi = 1.1e3)        
        sampleDict.add("t1_3_points", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim/sms_3_points.root"]', lumi = 1.1e3)

        return [sampleDict,samples.susy17,samples.top17]

    def listOfSamples(self,params) :
        from supy.samples import specify
        return (#supy.samples.specify(names = ["T2tt_8.job351"])+
                supy.samples.specify(names = ["T2tt_500_0"],   nEventsMax = 20000)+
                supy.samples.specify(names = ["T2tt_500_100"], nEventsMax = 20000)+
                supy.samples.specify(names = ["T2tt_500_300"], nEventsMax = 20000)+
                supy.samples.specify(names = ["tt_8_mg.job315_zeroNu_50Met"])+
                supy.samples.specify(names = ["tt_8_mg.job315_1_oneNu_100Met"])+
                supy.samples.specify(names = ["tt_8_mg.job315_1_twoNu_100Met"])+
                #supy.samples.specify(names = ["tt_8_mg.job315_1"])+
                []
                )

    def conclude(self,pars) :
        org = self.organizer(pars)

        def md(x, y) :
            x.update(y)
            return x
        mcOps = {"markerStyle":1, "lineWidth":1, "goptions":"ehist"}
        org.mergeSamples(targetSpec = md({"name":"st-st (500,  0)", "color":r.kGreen}, mcOps), allWithPrefix = "T2tt_500_0")
        org.mergeSamples(targetSpec = md({"name":"st-st (500,100)", "color":r.kRed}, mcOps), allWithPrefix = "T2tt_500_100")
        org.mergeSamples(targetSpec = md({"name":"st-st (500,300)", "color":r.kBlue}, mcOps), allWithPrefix = "T2tt_500_300")
        org.mergeSamples(targetSpec = md({"name":"tt", "color":r.kCyan}, mcOps), sources = ["tt_8_mg.job315_1"])
        org.mergeSamples(targetSpec = md({"name":"tt (0nu,50met)", "color":r.kBlack}, mcOps), sources = ["tt_8_mg.job315_zeroNu_50Met"])
        org.mergeSamples(targetSpec = md({"name":"tt (1nu,100met)", "color":r.kOrange+3}, mcOps), sources = ["tt_8_mg.job315_1_oneNu_100Met"])
        org.mergeSamples(targetSpec = md({"name":"tt (2nu,100met)", "color":r.kOrange-3}, mcOps), sources = ["tt_8_mg.job315_1_twoNu_100Met"])

        org.scale(20.0e3)
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          #doLog = False,
                          pegMinimum = 0.1,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          ).plotAll()
