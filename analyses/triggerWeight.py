import supy,steps,calculables,samples,math

class triggerWeight(supy.analysis) :
    def mutriggers(self) :             # L1 prescaling evidence
        ptv = { 12   :(1,2,3,4,5),     # 7,8,11,12
                15   :(2,3,4,5,6,8,9), # 12,13
                20   :(1,2,3,4,5,7,8),
                24   :(1,2,3,4,5,7,8,11,12),
                24.21:(1,),
                30   :(1,2,3,4,5,7,8,11,12),
                30.21:(1,),
                40   :(1,2,3,5,6,9,10),
                40.21:(1,4,5),
                }
        return sum([[("HLT_Mu%d%s_v%d"%(int(pt),"_eta2p1" if type(pt)!=int else "",v),int(pt)+1) for v in vs] for pt,vs in sorted(ptv.iteritems())],[])

    def unreliableTriggers(self) :
        '''Evidence of L1 prescaling at these ostensible prescale values'''
        return { "HLT_Mu15_v9":(25,65),
                 "HLT_Mu20_v8":(30,),
                 "HLT_Mu24_v8":(20,25),
                 "HLT_Mu24_v11":(35,),
                 "HLT_Mu24_v12":(35,),
                 "HLT_Mu30_v8":(4,10),
                 "HLT_Mu30_v11":(4,20),
                 "HLT_Mu30_v12":(8,20),
                 "HLT_Mu40_v6":(4,),
                 "HLT_Mu40_v9":(10,),
                 "HLT_Mu40_v10":(10,)
                 }

    def parameters(self) :
        return {"muon" : ("muon","PF"),
                "triggers": self.mutriggers()
                }
    
    def listOfCalculables(self,pars) :
        return (supy.calculables.zeroArgs() +
                supy.calculables.fromCollections(calculables.muon,[pars["muon"]]) +
                [calculables.muon.Indices( pars["muon"], ptMin = 10, combinedRelIsoMax = 0.15),
                 calculables.muon.IndicesTriggering( pars['muon'] ),
                 calculables.vertex.ID(),
                 calculables.vertex.Indices(),
                 supy.calculables.other.abbreviation( "nVertexRatio", "nvr" ),
                 supy.calculables.other.abbreviation('muonTriggerWeightPF','tw')
                 ])
    
    def listOfSteps(self,pars) :
        return [
            supy.steps.printer.progressPrinter(),
            supy.steps.other.histogrammer("genpthat",200,0,1000,title=";#hat{p_{T}} (GeV);events / bin"),
            supy.steps.filters.multiplicity("vertexIndices",min=1),
            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0"),
            steps.trigger.physicsDeclaredFilter(),
            steps.filters.monster(),
            steps.filters.hbheNoise(),
            calculables.trigger.TriggerWeight(samples = [ss.weightedName for ss in self.listOfSamples(pars) if 'SingleMu' in ss.name],
                                              triggers = zip(*pars['triggers'])[0], thresholds = zip(*pars['triggers'])[1],
                                              unreliable = self.unreliableTriggers()),
            supy.calculables.other.Ratio("nVertex", binning = (15,-0.5,14.5), thisSample = pars['baseSample'],
                                         target = ("SingleMu",[]), groups = [('qcd_py6',[]),('w_jets_fj_mg',[]),('tt_tauola_fj_mg',[])]),
            steps.histos.value("%sCombinedRelativeIso%s"%pars['muon'], 100, 0, 1, "%sIndicesTriggering%s"%pars['muon'], index=0),
            steps.histos.absEta("%sP4%s"%pars['muon'], 200,0,2.2, "%sIndicesTriggering%s"%pars['muon'], index=0),
            steps.histos.phi("%sP4%s"%pars['muon'], 200,-math.pi,math.pi, "%sIndicesTriggering%s"%pars['muon'], index=0),
            steps.histos.generic(("%sP4%s"%pars['muon'],"%sTriggeringIndex%s"%pars['muon']),
                                 (100,100),(0,-math.pi),(2.2,math.pi), title = ";mu |eta|;mu phi;events / bin",
                                 funcString = "lambda x: (abs(x[0][x[1]].eta()),x[0][x[1]].phi())" ),

            supy.steps.filters.absEta( "%sP4%s"%pars['muon'], min = 1.8, indices = "%sIndicesTriggering%s"%pars['muon'], index = 0),
            #supy.steps.filters.value("%sTriggeringPt%s"%pars['muon'], min = 35, max = 45),

            steps.histos.pt("%sP4%s"%pars['muon'], 200,0,200, "%sIndicesTriggering%s"%pars['muon'], index=0),
            steps.histos.pt("%sP4%s"%pars['muon'], 200,30,50, "%sIndicesTriggering%s"%pars['muon'], index=0),

            steps.histos.value("%sCombinedRelativeIso%s"%pars['muon'], 100, 0, 1, "%sIndicesTriggering%s"%pars['muon'], index=0),
            steps.histos.absEta("%sP4%s"%pars['muon'], 200,0,2.2, "%sIndicesTriggering%s"%pars['muon'], index=0),
            steps.histos.phi("%sP4%s"%pars['muon'], 200,-math.pi,math.pi, "%sIndicesTriggering%s"%pars['muon'], index=0),
            steps.histos.generic(("%sP4%s"%pars['muon'],"%sTriggeringIndex%s"%pars['muon']),
                                 (100,100),(0,-math.pi),(2.2,math.pi), title = ";mu |eta|;mu phi;events / bin",
                                 funcString = "lambda x: (abs(x[0][x[1]].eta()),x[0][x[1]].phi())" ),
            ]
            
    def listOfSampleDictionaries(self) : return [samples.muon,samples.mc]

    def listOfSamples(self,params) :
        return ( supy.samples.specify( names = ['SingleMu.2011B-PR1.1b',
                                                'SingleMu.2011B-PR1.1a',
                                                'SingleMu.2011A-Oct.1',
                                                'SingleMu.2011A-Aug.1',
                                                'SingleMu.2011A-PR4.1',
                                                'SingleMu.2011A-May.1'], weights = 'tw') +
                 supy.samples.specify( names = ['qcd_py6fjmu_pt_15_20',
                                                'qcd_py6fjmu_pt_20_30',
                                                'qcd_py6fjmu_pt_30_50',
                                                'qcd_py6fjmu_pt_50_80',
                                                'qcd_py6fjmu_pt_80_120',
                                                'qcd_py6fjmu_pt_120_150',
                                                'qcd_py6fjmu_pt_150',
                                                'w_jets_fj_mg',
                                                'tt_tauola_fj_mg'
                                                ], weights = ['tw','nvr'] , effectiveLumi = 300))
    
    def conclude(self,pars) :
        org = self.organizer(pars)
        org.mergeSamples(targetSpec = {"name":"SingleMu","color":r.kBlack,"markerStyle":20}, allWithPrefix="SingleMu")
        org.mergeSamples(targetSpec = {"name":"qcd", "color":r.kBlue}, allWithPrefix="qcd")
        org.mergeSamples(targetSpec = {"name":"w_jets","color":r.kRed}, allWithPrefix="w_jets")
        org.mergeSamples(targetSpec = {"name":"t#bar{t}","color":r.kViolet}, allWithPrefix="tt_tauola_fj_mg")
        #org.scaleOneRaw([ss['name'] for ss in org.samples].index('w_jets'), 0.6)
        org.mergeSamples(targetSpec = {"name":"s.m.", "color":r.kGreen+3}, keepSources = True, sources = ['qcd','w_jets','t#bar{t}'], force = True)
        org.scale()

        kwargs = { "blackList":["lumiHisto","xsHisto","nJobsHisto","muonTriggerWeightPF"],
                   "samplesForRatios":("SingleMu","s.m.") if "s.m." in [ss['name'] for ss in org.samples] else ("","")}

        supy.plotter(org, psFileName = self.psFileName(org.tag+'_log'),   doLog = True,  **kwargs).plotAll()
        supy.plotter(org, psFileName = self.psFileName(org.tag+'_nolog'), doLog = False, **kwargs).plotAll()
