import supy,steps,calculables,samples

class topAsymmTemplates(supy.analysis) :
    def parameters(self) :
        #weightsQQ = [calculables.top.wQQbarHardAsym(target = 0.0, nominal = a*0.05) for a in range(-20,21)]
        weightsQQ = [calculables.top.wQQbarHardFlat(a*0.05) for a in range(-20,21)]
        weightsQQName = [w.name for w in weightsQQ]
        #weightsQg = [calculables.top.wQgHardAsym(target = 0.0, nominal = a*0.05) for a in range(-20,21)]
        weightsQg = [calculables.top.wQgHardFlat(a*0.15) for a in range(-20,21)]
        weightsQgName = [w.name for w in weightsQg]
        return {"effectiveLumi" : 100000,
                "generator" : self.vary({"compare":["_mg","_ph","_mn"],
                                         "mg":"_mg",
                                         "ph":"_ph",
                                         "mn":"_mn"
                                         }),
                "weightsQQ" : weightsQQ,
                "weightsQQName" : weightsQQName,
                'weightsQg' : weightsQg,
                'weightsQgName' : weightsQgName
                }

    def listOfCalculables(self, pars) :
        return ( sum([supy.calculables.zeroArgs(module) for module in [calculables,supy.calculables]],[]) +
                 supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")]) +
                 [ calculables.vertex.ID(),
                   calculables.vertex.Indices(),
                   calculables.gen.genIndicesHardPartons( {"ttj_mg":(4,5),
                                                           "ttj_ph":(4,5),
                                                           "ttj_mn":(0,1)}[pars["baseSample"]] )
                   ]+pars["weightsQQ"]+pars["weightsQg"]
                 )
    
    def listOfSteps(self, pars) :
        return [supy.steps.printer.progressPrinter(),
                supy.calculables.other.SymmAnti(pars['sample'],"genttCosThetaStar",1),
                #steps.gen.topPrinter(),
                #supy.steps.histos.weighted("genttCosThetaStar", 50,-1,1, weights = pars["weightsQQName"], pred = "wQQ"),
                #supy.steps.histos.weighted("genttCosThetaStarBar", 50,-1,1, weights = pars["weightsQQName"], pred = "wQQ"),
                #supy.steps.histos.weighted("genttCosThetaStar", 100,-1,1, weights = pars["weightsQgName"], pred = "wQG"),
                #supy.steps.histos.weighted("genttCosThetaStarBar", 100,-1,1, weights = pars["weightsQgName"], pred = "wQG"),
                #supy.steps.filters.label("end templates"),
                #steps.gen.particlePrinter(),
                #steps.top.collisionType(),
                supy.steps.filters.value('wGG', max=0, allowNone = True),
                #steps.top.mcQuestions(),
                steps.top.mcQuestions2(pred = 'wQQ'),
                steps.top.mcQuestions2(pred = 'wQG'),
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
                                             color = col) for suf,col in suffixColor for w in ['wQQ','wQG']],[])

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
        org.scale(toPdf=True)

        from supy import plotter
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          doLog = False,
                          pegMinimum = 4e-4,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          ).plotAll()
