import supy,steps,calculables,samples
import os,math,copy,ROOT as r, numpy as np

class muIso(supy.analysis) :
    def parameters(self) :

        objects = {
            'label'       :[               'pf' ,               'pat' ],
            'jet'         :[  ('ak5JetPF','Pat'),     ('ak5Jet','Pat')],
            'muon'        :[       ('muon','PF'),       ('muon','Pat')],
            }

        return { "vary" : ['objects'],
                 "objects": self.vary([ ( objects['label'][index], dict((key,val[index]) for key,val in objects.iteritems())) for index in range(2) if objects['label'][index] in ['pf']]),
                 }

    def listOfSampleDictionaries(self) : return [getattr(samples,item) for item in ['muon16', 'top16', 'ewk16', 'qcd16']]

    def data(self,pars) :
        return supy.samples.specify( names = ['SingleMu.2011A.1',
                                              'SingleMu.2011A.2',
                                              'SingleMu.2011B'], weights = 'tw')

    def listOfSamples(self,pars) :

        def ewk(eL = None) :
            return supy.samples.specify( names = ["wj_lv_mg","dyj_ll_mg"], effectiveLumi = eL, color = 28 , nFilesMax = 20 )

        def ttbar_mg(eL = None) :
            return supy.samples.specify( names = "ttj_mg", effectiveLumi = eL, color = r.kBlue , nFilesMax = 20 )

        def qcd_py6_mu(eL = None) :
            return supy.samples.specify( names = ["qcd_mu_15_20",
                                                  "qcd_mu_20_30",
                                                  "qcd_mu_30_50",
                                                  "qcd_mu_50_80",
                                                  "qcd_mu_80_120",
                                                  "qcd_mu_120_150",
                                                  "qcd_mu_150"], nFilesMax = 20 )
        
        return  ( ewk() + ttbar_mg(5e4) + qcd_py6_mu() )

    ########################################################################################
    def listOfCalculables(self, pars) :
        obj = pars["objects"]
        calcs  = supy.calculables.zeroArgs(supy.calculables)
        calcs += supy.calculables.zeroArgs(calculables)
        for item in ['jet','muon'] :
            calcs += supy.calculables.fromCollections( getattr(calculables, item), [obj[item]])
        calcs += [
            calculables.jet.Indices(       obj["jet"],      ptMin = 20, etaMax = 3.0, flagName = "JetIDloose"),
            calculables.muon.Indices(      obj["muon"],     ptMin = 10, combinedRelIsoMax = 0.15),
            calculables.muon.IndicesAnyIsoIsoOrder(obj["muon"], "CombinedRelativeIso"),

            calculables.vertex.ID(),
            calculables.vertex.Indices(),
            ]
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :
        obj = pars["objects"]
        
        ssteps = supy.steps
        
        return (
            [ssteps.printer.progressPrinter()
             , ssteps.histos.value("genpthat",200,0,1000,xtitle="#hat{p_{T}} (GeV)").onlySim()
             , ssteps.histos.value("genQ",200,0,1000,xtitle="genQ (GeV)").onlySim()
             
             ####################################
             , ssteps.filters.label('data cleanup'),
             ssteps.filters.multiplicity("vertexIndices",min=1)

             , steps.other.muIsoStudy(obj['jet'],obj['muon'])

             ])

    ########################################################################################
    def concludeAll(self) :
        self.rowcolors = 2*[13] + 2*[45]
        super(muIso,self).concludeAll()

    def conclude(self,pars) :
        org = self.organizer(pars, verbose = True )
        org.mergeSamples(targetSpec={"name":"qcd"}, allWithPrefix="qcd")
        org.scale( toPdf = True )

        names = [ss["name"] for ss in org.samples]
        kwargs = {"detailedCalculables": False,
                  "blackList":["lumiHisto","xsHisto","nJobsHisto"],
                  "detailedCalculables" : True,
                  "rowColors" : self.rowcolors,
                  "rowCycle" : 100,
                  }
        
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_log"),  doLog = True, pegMinimum = 0.01, **kwargs ).plotAll()
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_nolog"), doLog = False, **kwargs ).plotAll()

