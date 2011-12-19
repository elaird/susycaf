import supy,steps,calculables,samples
import os,math,copy,ROOT as r, numpy as np

class kinfit(supy.analysis) :

    def parameters(self) :

        objects = {
            'label'       :[               'pf' ,               'pat' ],
            'muonsInJets' :[               True ,               False ],
            'jet'         :[('ak5JetPF','Pat'),   ('xcak5Jet','Pat')],
            'muon'        :[       ('muon','PF'),       ('muon','Pat')],
            'electron'    :[   ('electron','PF'),   ('electron','Pat')],
            'photon'      :[    ('photon','Pat'),     ('photon','Pat')],
            'met'         :[          'metP4PF' ,    'metP4AK5TypeII' ],
            'sumP4'       :[          'pfSumP4' ,           'xcSumP4' ],
            'sumPt'       :[       'metSumEtPF' ,           'xcSumPt' ],
            }

        return { "bVar" : "NTrkHiEff", # "TrkCountingHighEffBJetTags"
                 "objects": dict((key,val[0]) for key,val in objects.items()),
                 "lepton" : {'name':'muon'},
                 "topBsamples": ("ttj_mg",[]),
                 "lumi" : None,
                 "nFilesMax" : self.vary({'test':1, 'full':200})
                 }

    ########################################################################################

    def listOfSampleDictionaries(self) : return [samples.top]

    def listOfSamples(self,pars) :  return supy.samples.specify( names = pars['topBsamples'][0], effectiveLumi = pars['lumi'] , nFilesMax = pars['nFilesMax'], color = r.kRed )

    ########################################################################################
    def listOfCalculables(self, pars) :
        obj = pars["objects"]
        lepton = obj[pars["lepton"]["name"]]
        calcs  = supy.calculables.zeroArgs(supy.calculables)
        calcs += supy.calculables.zeroArgs(calculables)
        for item in ['jet','photon','electron','muon'] :
            calcs += supy.calculables.fromCollections( getattr(calculables, item), [obj[item]])
        calcs += [
            calculables.jet.IndicesBtagged(obj["jet"],pars["bVar"]),
            calculables.jet.Indices(       obj["jet"],      ptMin = 20, etaMax = 3.5, flagName = "JetIDloose"),
            calculables.photon.Indices(    obj["photon"],   ptMin = 25, flagName = "photonIDLooseFromTwikiPat"),
            calculables.electron.Indices(  obj["electron"], ptMin = 10, simpleEleID = "80", useCombinedIso = True),
            calculables.muon.Indices(      obj["muon"],     ptMin = 10, combinedRelIsoMax = 0.15),

            calculables.xclean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.xcJet(obj["jet"], applyResidualCorrectionsToData = False,
                                     gamma    = obj["photon"],      gammaDR = 0.5,
                                     electron = obj["electron"], electronDR = 0.5,
                                     muon     = obj["muon"],         muonDR = 0.5, correctForMuons = not obj["muonsInJets"]),

            calculables.top.mixedSumP4(transverse = obj["met"], longitudinal = obj["sumP4"]),
            calculables.top.SemileptonicTopIndex(lepton),
            calculables.top.fitTopLeptonCharge(lepton),
            calculables.top.TopReconstruction(lepton,obj["jet"],"mixedSumP4"),
            calculables.top.IndicesGenTopPQHL(obj["jet"]),
            calculables.top.genTopRecoIndex(jets=obj['jet']),
            calculables.top.genTopSemiLeptonicWithinAcceptance(jetPtMin = 20, jetAbsEtaMax = 3.5, lepPtMin = 20, lepAbsEtaMax = 2.1),

            calculables.top.TopComboQQBBLikelihood(pars['objects']['jet'], pars['bVar']),
            calculables.top.OtherJetsLikelihood(pars['objects']['jet'], pars['bVar']),

            calculables.other.Covariance(('met','PF')),

            supy.calculables.other.abbreviation( "TrkCountingHighEffBJetTags", "NTrkHiEff", fixes = calculables.jet.xcStrip(obj['jet']) ),
            ]
        calcs += supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")])
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :
        
        obj = pars["objects"]

        el = obj['electron']
        mu = obj['muon']
        jet = obj['jet']
        
        bVar = ("%s"+pars["bVar"]+"%s")%calculables.jet.xcStrip(jet)
        
        ssteps = supy.steps
        
        return (
            [ssteps.printer.progressPrinter(),
             ssteps.filters.value("genTopSemiMu", min=1),
             ssteps.filters.value("genTopSemiLeptonicWithinAcceptance", min=1)
             
             , ssteps.filters.label('selection'),

             ssteps.filters.pt(obj['met'],min=20),
             ssteps.filters.multiplicity("%sIndices%s"%el, max=0       ),
             ssteps.filters.multiplicity("%sIndices%s"%mu, max=1, min=1),
             ssteps.filters.absEta("%sP4%s"%mu, max = 2.1,  indices = "%sIndices%s"%mu, index = 0),
             ssteps.filters.pt(    "%sP4%s"%mu, min = 20.0, indices = "%sIndices%s"%mu, index = 0),
             ssteps.filters.multiplicity("%sIndices%s"%jet, min = 4 )

             , calculables.jet.ProbabilityGivenBQN(jet, pars['bVar'], binning=(64,-1,15), samples = pars['topBsamples']),
             ssteps.filters.value(bVar, indices = "%sIndicesBtagged%s"%jet, index=1, min=2.0)

             , ssteps.filters.label('top reco'),
             ssteps.filters.multiplicity("TopReconstruction",min=1),

             ssteps.filters.minimum("%sIndicesGenTopPQHL%s"%jet, 0),
             ssteps.filters.unique("%sIndicesGenTopPQHL%s"%jet)

             #, steps.top.jetPrinter(jet)
             
             , ssteps.filters.label("fitTop kinFitLook")            , steps.top.kinFitLook("fitTopRecoIndex")
             , ssteps.filters.label("genTop kinFitLook")            , steps.top.kinFitLook("genTopRecoIndex")
             , ssteps.filters.label("combinatorial filtering")      , steps.top.combinatorialFiltering(jet)
             , ssteps.filters.label("combinatorial frequency")      , steps.top.combinatorialFrequency(jet)
             , ssteps.filters.label("genTopRecoIndex resolutions")  , steps.top.resolutions('genTopRecoIndex')
             , ssteps.filters.label("fitTopRecoIndex resolutions")  , steps.top.resolutions('fitTopRecoIndex')
             ])

    ########################################################################################

    def conclude(self,pars) :
        org = self.organizer(pars, verbose = True )
        org.scale( toPdf = True )
        #org.scale( lumiToUseInAbsenceOfData = pars['lumi'] )

        kwargs = {"detailedCalculables": False,
                  "blackList":["lumiHisto","xsHisto","nJobsHisto"],
                  "rowColors" : 2*[13] + 2*[45],
                  "rowCycle" : 100,
                  #"omit2D" : True,
                  }
        
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_log"),  doLog = True, pegMinimum = 0.01, **kwargs ).plotAll()
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_nolog"), doLog = False, **kwargs ).plotAll()
        r.gStyle.SetOptStat(0)

        self.resPlots(org)
        self.filterPlots(org)
        self.fitPlots(org)

    def resPlots(self,org) :
        genRes = org.steps[next(org.indicesOfStep('resolutions','genTopRecoIndex'))]
        fitRes = org.steps[next(org.indicesOfStep('resolutions','fitTopRecoIndex'))]
        
        for var in  ["Rapidity",'eta'] :
            for subvar,title in [(s%var,t%var) for s,t in [('d%sLepTop_','#Delta %s_{reco-gen} leptonic top'),
                                                           ('d%sHadTop_','#Delta %s_{reco-gen} hadronic top'),
                                                           ('dd%sTTbar_','#Delta_{reco-gen} #Delta %s_{t#bar{t}}')]] :
                fit = subvar+'fit'
                raw = subvar+'unfit'
                hists = [
                    ("trueComb fit" , genRes[fit][0], r.kBlack, r.kSolid, 2),
                    ("trueComb raw" , genRes[raw][0], r.kBlack, r.kDashed, 1),
                    ("selected fit" , fitRes[fit][0], r.kRed, r.kSolid, 2),
                    ("selected raw" , fitRes[raw][0], r.kRed, r.kDashed, 1)
                    ]

                m = max( hist.GetMaximum() for n,hist,c,d,w in hists)
                L = r.TLegend(0.1,0.7,0.38,0.9)
                can = r.TCanvas()
                for i,(n,hist,c,d,w) in enumerate(hists) :
                    hist.SetMaximum(1.1*m)
                    hist.SetLineWidth(w)
                    hist.SetLineColor(c)
                    hist.SetLineStyle(d)
                    hist.SetTitle(";%s;"%title)
                    L.AddEntry(hist,n,"l")
                    hist.Draw('histsame' if i else 'hist')
                L.Draw()
                can.Update()
                can.Print(subvar+'.pdf', 'pdf')

        
    def filterPlots(self,org) :
        fltr = org.steps[next(org.indicesOfStep('combinatorialFiltering'))]
        can = r.TCanvas()

        fltr['max_bjet_genIndex'][0].SetLineWidth(2)
        fltr['max_bjet_genIndex'][0].SetLineColor(r.kBlue)
        fltr['max_bjet_genIndex'][0].SetTitle(';reco b-index of gen b with lesser TrackCountingHE;')
        fltr['max_bjet_genIndex'][0].Draw('hist')
        cut1 = r.TLine(4.5,0,4.5,0.7)
        cut2 = r.TArrow(4.5,0.7,3.5,0.7)
        cut1.SetLineWidth(2)
        cut2.SetLineWidth(2)
        cut1.Draw()
        cut2.Draw()
        can.Print('max_bjet_genIndex.pdf','pdf')


        can.Clear()
        can.Divide(2,2)
        for i,(c,p) in enumerate([(c,p) for p in ['','pass'] for c in ['correct','incorrect']]) :
            print i,c,p
            can.cd(i+1)
            fltr["combo_raw_m_%s_%s"%(c,p)][0].SetTitle('%s jet comb.'%c+(" passing cut" if p else ''))
            fltr["combo_raw_m_%s_%s"%(c,p)][0].Draw('colz')
        can.Print("combo_raw_m.pdf","pdf")


    def fitPlots(self,org) :
        gen = org.steps[next(org.indicesOfStep("kinFitLook","genTopRecoIndex"))]
        fit = org.steps[next(org.indicesOfStep("kinFitLook","fitTopRecoIndex"))]

        for (f,fH),(g,gH) in zip(sorted(fit.items()),sorted(gen.items())) :
            can = r.TCanvas()
            fH[0].SetLineColor(r.kRed)
            fH[0].SetLineWidth(2)
            gH[0].SetLineWidth(2)
            gH[0].SetLineColor(r.kBlack)
            fH[0].SetTitle(';%s;'%fH[0].GetXaxis().GetTitle().replace("fitTopRecoIndex",""))
            gH[0].SetTitle(';%s;'%gH[0].GetXaxis().GetTitle().replace("genTopRecoIndex",""))
            fH[0].Draw("hist")
            gH[0].Draw("histsame")
            L = r.TLegend(0.58,0.7,0.9,0.9)
            L.AddEntry(fH[0],'selected',"l")
            L.AddEntry(gH[0],'true',"l")
            L.Draw()
            can.Print(f.replace("fitTopRecoIndex","")+".pdf",'pdf')
            
