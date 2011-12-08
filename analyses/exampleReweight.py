import supy,samples,calculables,ROOT as r

class exampleReweight(supy.analysis) :

    def listOfSteps(self,pars) :
        return [ supy.steps.printer.progressPrinter(),
                 supy.steps.histos.value("genpthat",200,0,1000,xtitle=";#hat{p_{T}} (GeV)").onlySim(),
                 supy.steps.histos.multiplicity("vertexIndices",max=15),
                 supy.steps.filters.multiplicity("vertexIndices",min=1),
                 supy.steps.filters.pt("muonP4PF", min = 25, indices = "muonIndicesPF", index=0),
                 supy.steps.histos.multiplicity("vertexIndices",max=15),
                 supy.calculables.other.Ratio("nVertex", binning = (15,-0.5,14.5), thisSample = pars['baseSample'],
                                              target = ("SingleMu",[]), groups = [('qcd_mg',[]),('qcd_py6',[])],
                                              ),
                 supy.steps.histos.multiplicity("vertexIndices",max=15),
                 ]
    
    def listOfCalculables(self,pars) :
        muon = ("muon","PF")
        return ( supy.calculables.zeroArgs(supy.calculables) +
                 supy.calculables.zeroArgs(calculables) +
                 supy.calculables.fromCollections(calculables.muon, [muon]) +
                 [ calculables.muon.Indices( muon, ptMin = 10, combinedRelIsoMax = 0.15),
                   calculables.vertex.ID(),
                   calculables.vertex.Indices(),
                   ] )

    def listOfSampleDictionaries(self) :
        return [samples.mc,samples.muon]

    def listOfSamples(self,pars) :
        def qcd_mg(eL) :
            qM = ["%d"%t for t in [50,100,250,500,1000][1:]]
            return supy.samples.specify( effectiveLumi = eL, weights = "nVertexRatio",
                                         names = ["qcd_mg_ht_%s_%s"%t for t in zip(qM,qM[1:]+["inf"])])
        def qcd_py6_mu(eL) :
            q6 = [0,5,15,20,30,50,80,120,150,None]
            iCut = q6.index(15)
            return supy.samples.specify( effectiveLumi = eL, weights = "nVertexRatio",
                                         names = ["qcd_py6fjmu_pt_%s"%("%d_%d"%(low,high) if high else "%d"%low)
                                                  for low,high in zip(q6[:-1],q6[1:])[iCut:]] )
        return (supy.samples.specify(names="SingleMu.Run2011A-PR-v4.FJ.Burt4") +
                qcd_mg(100) +
                qcd_py6_mu(100) +
                [])
                
    def conclude(self,pars) :
        org = self.organizer(pars)
        org.mergeSamples(targetSpec = {"name":"qcd_mg", "color":r.kBlue}, allWithPrefix="qcd_mg")
        org.mergeSamples(targetSpec = {"name":"qcd_py6", "color":r.kRed}, allWithPrefix="qcd_py6")
        org.scale()
        supy.plotter( org,
                      pdfFileName = self.pdfFileName(org.tag),
                      blackList = ["lumiHisto","xsHisto","xsPostWeightsHisto","nJobsHisto","genpthat"],
                      detailedCalculables = True,
                      ).plotAll()
        
