import supy,steps,calculables,samples,ROOT as r

class topAsymmTemplates(supy.analysis) :
    def parameters(self) :
        #weightsQQ = [calculables.top.wQQbarHardAsym(target = 0.0, nominal = a*0.05) for a in range(-20,21)]
        #weightsQQ = [calculables.top.wQQbarHardFlat(a*0.05) for a in range(-20,21)]
        #weightsQQName = [w.name for w in weightsQQ]
        #weightsQg = [calculables.top.wQgHardAsym(target = 0.0, nominal = a*0.05) for a in range(-20,21)]
        #weightsQg = [calculables.top.wQgHardFlat(a*0.15) for a in range(-20,21)]
        #weightsQgName = [w.name for w in weightsQg]
        return {"effectiveLumi" : 100000,
                "generator" : self.vary({"compare":["_mg","_ph","_mn"],
                                         "mg":"_mg",
                                         "ph":"_ph",
                                         "mn":"_mn"
                                         }),
                #"weightsQQ" : weightsQQ,
                #"weightsQQName" : weightsQQName,
                #'weightsQg' : weightsQg,
                #'weightsQgName' : weightsQgName
                }

    def listOfCalculables(self, pars) :
        return ( sum([supy.calculables.zeroArgs(module) for module in [calculables,supy.calculables]],[]) +
                 supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")]) +
                 [ calculables.vertex.ID(),
                   calculables.vertex.Indices(),
                   calculables.gen.genIndicesHardPartons( {"ttj_mg":(4,5),
                                                           "ttj_ph":(4,5),
                                                           "ttj_mn":(0,1)}[pars["baseSample"]] )
                   ]#+pars["weightsQQ"]+pars["weightsQg"]
                 )
    
    def listOfSteps(self, pars) :
        return [supy.steps.printer.progressPrinter(),
                #supy.calculables.other.SymmAnti(pars['sample'],"genTopCosThetaBoost",1, inspect=True),
                supy.calculables.other.SymmAnti(pars['sample'],"genTopCosPhiBoost",1, inspect=True, nbins=160,
                                                funcEven = r.TF1('phiboost',"[0]*(1+[1]*x**2)/sqrt(1-x**2)",-1,1),
                                                funcOdd = r.TF1('phiboostodd','[0]*x/sqrt(1-x**2)',-1,1)),
                supy.calculables.other.SymmAnti(pars['sample'],"genttCosThetaStar",1, inspect=True,
                                                funcEven = '++'.join('x**%d'%(2*d) for d in range(5)),
                                                funcOdd = '++'.join('x**%d'%(2*d+1) for d in range(5))),
                supy.calculables.other.SymmAnti(pars['sample'],"genTopCosThetaBoostAlt",1, inspect=True,
                                                funcEven = '++'.join('x**%d'%(2*d) for d in range(5)),
                                                funcOdd = '++'.join('x**%d'%(2*d+1) for d in range(5))),
                supy.calculables.other.SymmAnti(pars['sample'],"genTopDeltaBetazRel",1, inspect=True,
                                                funcEven = '++'.join(['(1-abs(x))']+['x**%d'%d for d in [0,2,4,6,8,10,12,14,16,18]]),
                                                funcOdd = '++'.join(['x**%d'%d for d in [1,3,5,7,9,11,13]])),
                #supy.calculables.other.SymmAnti(pars['sample'],"genTopDeltaAbsYttbar",3, inspect=True),
                #supy.steps.filters.label('reweighting'),
                #supy.steps.histos.symmAnti('hard','genttOm',100,-2,2),
                #supy.steps.histos.symmAnti('hard','genCosThetaStar',100,-1,1),
                #supy.steps.histos.symmAnti('hard','genCosThetaStarBar',100,-1,1),
                #supy.steps.histos.symmAnti('hard','genttCosThetaStar',100,-1,1),
                #supy.steps.histos.symmAnti('hard','genttCosThetaStarBar',100,-1,1),
                #supy.steps.histos.symmAnti('hard','genTopCosThetaBoost',100,-1,1),
                #steps.gen.topPrinter(),
                #supy.steps.histos.weighted("genttCosThetaStar", 50,-1,1, weights = pars["weightsQQName"], pred = "wQQ"),
                #supy.steps.histos.weighted("genttCosThetaStarBar", 50,-1,1, weights = pars["weightsQQName"], pred = "wQQ"),
                #supy.steps.histos.weighted("genttCosThetaStar", 100,-1,1, weights = pars["weightsQgName"], pred = "wQG"),
                #supy.steps.histos.weighted("genttCosThetaStarBar", 100,-1,1, weights = pars["weightsQgName"], pred = "wQG"),
                #supy.steps.filters.label("end templates"),
                #steps.gen.particlePrinter(),
                #steps.top.collisionType(),
                #supy.steps.filters.value('wGG', max=0, allowNone = True),
                #steps.top.mcQuestions(),
                steps.top.mcQuestions2(),
                #steps.filters.label("all"),         steps.Top.mcTruthTemplates(),
                #steps.filters.OR([steps.Filter.value('genTTbarIndices',min=0,index='lplus'),
                #                 steps.Filter.value('genTTbarIndices',min=0,index='lminus')]),
                #steps.top.mcTruthTemplates(),
                #steps.filters.label("acceptance"),        steps.Top.mcTruthAcceptance(),
                #steps.filters.label("discriminateQQbar"), steps.Top.discriminateQQbar(('genTop','')),
                #steps.filters.label("q direction"),       steps.Top.mcTruthQDir(),
                ]
    
    def listOfSampleDictionaries(self) : return [samples.top16]

    def listOfSamples(self,pars) :
        import ROOT as r
        eL = pars["effectiveLumi"]

        if type(pars["generator"]) is list :
            suffixColor = zip(pars["generator"],[r.kBlack,r.kRed,r.kBlue])
            return sum([supy.samples.specify(names = "ttj%s"%suf, effectiveLumi = eL, weights = [w],
                                             color = col + (0 if w=='wQQ' else 2)) for suf,col in suffixColor for w in ['wQQ','wQG']],[])

        sample = "ttj%s"%pars["generator"]
        asymms = [(r.kBlue, -0.3),
                  (r.kGreen, 0.0),
                  (r.kRed,   0.3)]
        R_sm = -0.05 if pars['generator'] == "mg" else 0.0
        return (
            #supy.samples.specify( names = sample, effectiveLumi = 500, color = r.kBlack,     weights = calculables.Gen.wNonQQbar()) +
            #supy.samples.specify( names = sample, effectiveLumi = eL, color = r.kRed,       weights = calculables.Gen.wQQbar()) +
            sum([supy.samples.specify(names = sample, nFilesMax = 4, #effectiveLumi = eL,
                                      color = col, weights = calculables.top.wTopAsym(R,R_sm=R_sm)) for col,R in asymms],[]) +
            [])
    
    def conclude(self,pars) :
        org = self.organizer(pars)
        org.scale(1000,toPdf=False)

        from supy import plotter
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          doLog = False,
                          pegMinimum = 4e-4,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          ).plotAll()
