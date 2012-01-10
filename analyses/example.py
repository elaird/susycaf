import supy,calculables,steps,samples, ROOT as r

class example(supy.analysis) :
    def parameters(self) :
        return {"etRatherThanPt" : [False],
                "jets" : ("ak5JetPF","Pat"),
                "jetPtMin" : 10.0,
                }

    def listOfSteps(self,pars) :
        jets = pars["jets"]
        pt = pars["jetPtMin"]
        
        outList=[
            supy.steps.printer.progressPrinter(),
            steps.trigger.techBitFilter([0],True).onlyData(),
            steps.trigger.physicsDeclaredFilter().onlyData(),
            supy.steps.histos.multiplicity("vertexIndices"),
            supy.steps.filters.multiplicity("vertexIndices", min = 1),
            steps.filters.monster(),

            supy.steps.filters.pt("%sCorrectedP4%s"%jets, min = pt, indices = "%sIndices%s"%jets, index = 0),
            supy.steps.filters.pt("%sCorrectedP4%s"%jets, min = pt, indices = "%sIndices%s"%jets, index = 1),
            supy.steps.filters.multiplicity("%sIndicesOther%s"%jets, max = 0),
            supy.steps.filters.multiplicity("%sIndices%s"%jets, min = 2),
            
            steps.jet.singleJetHistogrammer(jets,1), 
            steps.jet.cleanJetHtMhtHistogrammer(jets, pars["etRatherThanPt"]),
            #steps.other.variableGreaterFilter(25.0,jets[0]+"SumPt"+jets[1]),
            
            steps.jet.alphaHistogrammer(jets, etRatherThanPt = pars["etRatherThanPt"]),
            #steps.other.skimmer(),
            ]
        return outList
    
    def listOfCalculables(self,pars) :
        jets = pars["jets"]
        pt = pars["jetPtMin"]
        listOfCalculables = supy.calculables.zeroArgs(supy.calculables)
        listOfCalculables = supy.calculables.zeroArgs(calculables)
        listOfCalculables += supy.calculables.fromCollections(calculables.jet,[jets])
        listOfCalculables += [
            calculables.jet.Indices( jets, ptMin = pt, etaMax = 3.0, flagName = "JetIDloose"),
            calculables.jet.SumP4( jets),
            calculables.jet.DeltaPhiStar( jets ),
            calculables.jet.AlphaT        ( jets, pars["etRatherThanPt"]),
            calculables.jet.DeltaPseudoJet( jets, pars["etRatherThanPt"]),
            calculables.vertex.ID(),
            calculables.vertex.Indices(),
            ]
        return listOfCalculables

    def listOfSampleDictionaries(self) :
        exampleDict = supy.samples.SampleHolder()
        exampleDict.add("Example_Skimmed_900_GeV_Data", '["/afs/cern.ch/user/e/elaird/public/susypvt/framework_take3/skimmed_900_GeV_Data.root"]', lumi = 1.0e-5 ) #/pb
        exampleDict.add("Example_Skimmed_900_GeV_MC", '["/afs/cern.ch/user/e/elaird/public/susypvt/framework_take3/skimmed_900_GeV_MC.root"]',       xs = 1.0e8 ) #pb
        return [exampleDict]

    def listOfSamples(self,pars) :
        return (supy.samples.specify(names = "Example_Skimmed_900_GeV_Data", color = r.kBlack, markerStyle = 20) +
                supy.samples.specify(names = "Example_Skimmed_900_GeV_MC", color = r.kRed, effectiveLumi = 0.5) )

    def conclude(self,pars) :
        #make a pdf file with plots from the histograms created above
        org = self.organizer(pars)
        org.scale()
        supy.plotter( org,
                      pdfFileName = self.pdfFileName(org.tag),
                      samplesForRatios = ("Example_Skimmed_900_GeV_Data","Example_Skimmed_900_GeV_MC"),
                      sampleLabelsForRatios = ("data","sim"),
                      ).plotAll()
