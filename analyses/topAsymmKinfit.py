import supy,topAsymmShell,calculables,samples

class topAsymmKinfit(topAsymmShell.topAsymmShell) :
    def parameters(self) :
        pars = super(topAsymmKinfit,self).parameters()
        pars["effectiveLumi"] = 8000
        pars["topBsamples"] = ("tt_tauola_fj_mg",["tt_tauola_fj_mg"])
        return pars

    def listOfCalculables(self,pars) :
        calcs = super(topAsymmKinfit,self).listOfCalculables(pars)
        calcs.append( calculables.top.genTopSemiLeptonicWithinAcceptance( jetPtMin = 20, jetAbsEtaMax=3.5, lepPtMin=31, lepAbsEtaMax = 2.1))
        calcs.append( calculables.top.genTopSemiLeptonicAccepted( pars['objects']['jet']))
        calcs.append( calculables.top.genTopRecoIndex())
        return calcs

    def listOfSteps(self, pars) :
        import steps
        obj = pars["objects"]
        lepton = obj[pars["lepton"]["name"]]
        lPtMin = pars["lepton"]["ptMin"]
        bVar = ("%s"+pars["bVar"]+"%s")%calculables.Jet.xcStrip(obj["jet"])
        
        return ([
            supy.steps.printer.progressPrinter(),
            supy.steps.filters.pt("%sP4%s"%lepton, min = lPtMin, indices = "%sIndicesAnyIso%s"%lepton, index = 0),
            ]+topAsymmShell.topAsymmShell.cleanupSteps(pars)+[
            ]+topAsymmShell.topAsymmShell.selectionSteps(pars, withPlots = False) +[
            #steps.top.kinFitLook("fitTopRecoIndex"),
            supy.steps.filters.value("genTopSemiLeptonicWithinAcceptance", min = True),
            #steps.Histos.value("genTopWqqDeltaR",50,0,4),
            supy.steps.filters.value("genTopSemiLeptonicAccepted", min = True),
            #steps.Histos.value("genTopWqqDeltaR",50,0,4),
            #steps.top.topProbLook(obj['jet']),
            supy.steps.other.assertNotYetCalculated("TopReconstruction"),
            supy.steps.filters.multiplicity("TopReconstruction",min=1),
            supy.steps.filters.value("genTopRecoIndex", min = 0),
            steps.top.combinatorialBG(obj['jet']),
            #steps.Histos.value("genTopWqqDeltaR",50,0,4),
            #steps.Histos.value("fitTopWqqDeltaR",50,0,4),
            supy.steps.filters.label('selected combo'), steps.top.kinFitLook("fitTopRecoIndex"), steps.top.combinatoricsLook("fitTopRecoIndex"),
            supy.steps.filters.label('true combo'),  steps.top.kinFitLook("genTopRecoIndex"), steps.top.combinatoricsLook("genTopRecoIndex", jets = obj['jet']),
            supy.steps.filters.multiplicity("%sIndices%s"%obj["jet"], min=4, max=4),
            steps.top.combinatorialBG(obj['jet']),
            steps.top.combinatoricsLook("genTopRecoIndex", jets = obj['jet']),
            steps.top.combinatoricsLook("fitTopRecoIndex"),
            ])
    
    def listOfSamples(self,pars) :
        import ROOT as r
        return (supy.samples.specify(names = "tt_tauola_fj_mg", color = r.kRed,
                                     #nFilesMax=1, nEventsMax=4100) +
                                     effectiveLumi = pars["effectiveLumi"]) +
                supy.samples.specify(names = "tt_tauola_fj_mg", color = r.kBlue, weights = "wQQbar") +
                #supy.samples.specify(names = "tt_tauola_fj", color = r.kBlue,
                #                     #nFilesMax=1, nEventsMax=4100)+
                #                       effectiveLumi = pars["effectiveLumi"]) +
                #supy.samples.specify( names = "w_jets_mg", effectiveLumi = 100, color = 28 ) +
                [])
    
    def conclude(self,pars) :
        org = self.organizer(pars)
        org.scale(toPdf=True)
        
        supy.plotter(org,
                     psFileName = self.psFileName(org.tag),
                     doLog = False,
                     #noSci = True,
                     #pegMinimum = 0.1,
                     detailedCalculables = True,
                     blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                     ).plotAll()

