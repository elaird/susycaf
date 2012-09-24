import supy,steps,calculables,samples
import os,math,copy,itertools,ROOT as r, numpy as np

class ttKinTest(supy.analysis) :
    '''Test kinematic fit of semileptonic top decay and kinematic solutions of dileptonic top decay using smeared gen-level objects.'''

    def parameters(self) :
        return { "vary" : ['toptype'],
                 "toptype" : self.vary({"mg":"mg","ph":"ph"}),
                 }

    def listOfSampleDictionaries(self) :
        return [samples.top16]

    def listOfSamples(self,pars) :
        return supy.samples.specify(names = "ttj_%s"%pars['toptype'], effectiveLumi = 10000, color=r.kBlue)

    ########################################################################################
    def listOfCalculables(self, pars) :
        calcs = sum( [supy.calculables.zeroArgs(module) for module in [calculables, supy.calculables]], [])
        calcs += supy.calculables.fromCollections(calculables.top,[('genTop',"")])
        calcs += [ calculables.gen.genIndicesHardPartons() ]
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :
        tt = pars['toptype']
        ssteps = supy.steps

        return (
            [ssteps.printer.progressPrinter(),
             ssteps.filters.value('genTTbarIndices', index = 'lplus', min = 0),
             ssteps.filters.value('genTTbarIndices', index = 'lminus', min = 0),
             ssteps.filters.label('smear'),
             steps.top.dileptonSolver(),
             ssteps.filters.label('gen'),
             steps.top.dileptonSolver(gen=True)
             ])
    ########################################################################################
    def concludeAll(self) :
        self.rowcolors = 2*[13] + 2*[45]
        super(ttKinTest,self).concludeAll()

    def conclude(self,pars) :
        org = self.organizer(pars, verbose = True )
        org.scale( lumiToUseInAbsenceOfData = 5008 )
        kwargs = {"detailedCalculables": True,
                  "blackList":["lumiHisto","xsHisto","nJobsHisto"],
                  "rowColors" : self.rowcolors,
                  "rowCycle" : 100,
                  }
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_log"),  doLog = True, pegMinimum = 0.01, **kwargs ).plotAll()
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_nolog"), doLog = False, **kwargs ).plotAll()
