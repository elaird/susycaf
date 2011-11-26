import supy,samples,calculables,steps

class topAsymmShell(supy.analysis) :

    def parameters(self) :
        def mutriggers() : # FIXME
            ptv = {
                #3:(3,4),
                #5:(3,4,5,6),
                #8:(1,2,3,4),
                12:(1,2,3,4,5),
                15:(2,3,4,5,6),
                20:(1,2,3,4,5),
                24:(1,2,3,4,5),
                30:(1,2,3,4,5),
                40:(1,2,3),
                #100:(1,2,3),
                }
            return sum([[("HLT_Mu%d_v%d"%(pt,v),pt+1) for v in vs] for pt,vs in sorted(ptv.iteritems())],[])
        
        objects = self.vary()
        fields =                           [ "jet",               "met",           "sumP4",    "sumPt",       "muon",         "electron",         "photon",         "muonsInJets"]
        objects["pf"]   = dict(zip(fields, [("xcak5JetPF","Pat"), "metP4PF",       "pfSumP4",  "metSumEtPF",  ("muon","PF"),  ("electron","PF"),  ("photon","Pat"),  True]))
        #objects["calo"] = dict(zip(fields, [("xcak5Jet","Pat"),  "metP4AK5TypeII", "xcSumP4", "xcSumPt",     ("muon","Pat"), ("electron","Pat"), ("photon","Pat"),  False]))

        leptons = self.vary()
        fieldsLepton    =                            ["name","ptMin", "etaMax",              "isoVar", "triggers"]
        leptons["muon"]     = dict(zip(fieldsLepton, ["muon",     20,     2.1, "CombinedRelativeIso",   mutriggers()]))
        #leptons["electron"] = dict(zip(fieldsLepton, ["electron", 30,       9,         "IsoCombined", ("FIX","ME")]))
        
        bVar = "NTrkHiEff" # "TrkCountingHighEffBJetTags"
        bCut = {"normal"   : {"index":1, "min":2.0},
                "inverted" : {"index":1, "max":2.0}}
        lIso = {"normal":  {"N":1, "indices":"Indices"},
                "inverted":{"N":0, "indices":"IndicesNonIso"}}
        
        return { "objects": objects,
                 "lepton" : leptons,
                 "nJets" :  {"min":4,"max":None},
                 "nJets2" : {"min":4,"max":None},
                 "bVar" : bVar,
                 "selection" : self.vary({"top" : {"bCut":bCut["normal"],  "lIso":lIso["normal"]},
                                          #"Wlv" : {"bCut":bCut["inverted"],"lIso":lIso["normal"]},
                                          "QCD" : {"bCut":bCut["normal"],  "lIso":lIso["inverted"]}
                                          })
                 }

    def listOfCalculables(self, pars) :
        obj = pars["objects"]
        lepton = obj[pars["lepton"]["name"]]
        outList  = supy.calculables.zeroArgs(supy.calculables)
        outList += supy.calculables.zeroArgs(calculables)
        outList += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])
        outList += supy.calculables.fromCollections(calculables.electron, [obj["electron"]])
        outList += supy.calculables.fromCollections(calculables.photon, [obj["photon"]])
        outList += supy.calculables.fromCollections(calculables.jet, [obj["jet"]])
        outList += [
            calculables.jet.IndicesBtagged(obj["jet"],pars["bVar"]),
            calculables.jet.Indices(      obj["jet"],      ptMin = 20, etaMax = 3.5, flagName = "JetIDloose"),
            calculables.muon.Indices(     obj["muon"],     ptMin = 10, combinedRelIsoMax = 0.15),
            calculables.muon.IndicesTriggering(obj["muon"]),
            calculables.electron.Indices( obj["electron"], ptMin = 10, simpleEleID = "80", useCombinedIso = True),
            calculables.photon.Indices(   obj["photon"],   ptMin = 25, flagName = "photonIDLooseFromTwikiPat"),

            calculables.xclean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.xcJet(obj["jet"], applyResidualCorrectionsToData = False,
                                     gamma    = obj["photon"],      gammaDR = 0.5,
                                     electron = obj["electron"], electronDR = 0.5,
                                     muon     = obj["muon"],         muonDR = 0.5, correctForMuons = not obj["muonsInJets"]),
            calculables.xclean.SumP4(obj["jet"], obj["photon"], obj["electron"], obj["muon"]),
            calculables.xclean.SumPt(obj["jet"], obj["photon"], obj["electron"], obj["muon"]),

            calculables.vertex.ID(),
            calculables.vertex.Indices(),
            calculables.other.lowestUnPrescaledTrigger(zip(*pars["lepton"]["triggers"])[0]),

            calculables.top.mixedSumP4(transverse = obj["met"], longitudinal = obj["sumP4"]),
            calculables.other.pt("mixedSumP4"),
            calculables.top.SemileptonicTopIndex(lepton),            
            calculables.top.fitTopLeptonCharge(lepton),
            calculables.top.TopReconstruction(lepton,obj["jet"],"mixedSumP4"),
            
            calculables.other.Mt(lepton,"mixedSumP4", allowNonIso=True, isSumP4=True),
            calculables.muon.IndicesAnyIsoIsoOrder(obj[pars["lepton"]["name"]], pars["lepton"]["isoVar"]),
            calculables.other.PtSorted(obj['muon']),
            calculables.other.Covariance(('met','PF')),
            supy.calculables.other.abbreviation( "TrkCountingHighEffBJetTags", "NTrkHiEff", fixes = calculables.jet.xcStrip(obj['jet']) ),
            supy.calculables.other.abbreviation( "nVertexRatio", "nvr" ),
            supy.calculables.other.abbreviation('muonTriggerWeightPF','tw'),
            calculables.jet.pt( obj['jet'], index = 0, Btagged = True ),
            calculables.jet.absEta( obj['jet'], index = 3, Btagged = False)
            ]
        outList += supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")])
        outList.append( calculables.top.TopComboQQBBLikelihood(pars['objects']['jet'], pars['bVar']))
        outList.append( calculables.top.OtherJetsLikelihood(pars['objects']['jet'], pars['bVar']))
        outList.append( calculables.top.TopRatherThanWProbability(priorTop=0.5) )
        return outList

    def listOfSampleDictionaries(self) :
        return [samples.mc, samples.muon]
    

    @staticmethod
    def dataCleanupSteps(pars) :
        obj = pars['objects']
        return ([
            steps.filters.hbheNoise(),
            steps.trigger.physicsDeclaredFilter(),
            steps.filters.monster(),
            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0"),
            steps.trigger.hltPrescaleHistogrammer(zip(*pars['lepton']['triggers'])[0]),
            steps.trigger.lowestUnPrescaledTriggerHistogrammer(),
            supy.steps.histos.multiplicity("vertexIndices", max=15),
            supy.steps.histos.value("%sPtSorted%s"%obj['muon'], 2,-0.5,1.5),
            supy.steps.filters.multiplicity("vertexIndices",min=1),
            ])

    @staticmethod
    def xcleanSteps(pars) :
        obj = pars['objects']
        return ([
            supy.steps.filters.multiplicity(s, max = 0) for s in ["%sIndices%s"%obj["photon"],
                                                                  "%sIndicesUnmatched%s"%obj["photon"],
                                                                  "%sIndices%s"%(obj["electron" if pars["lepton"]["name"]=="muon" else "muon"]),
                                                                  "%sIndicesUnmatched%s"%obj["electron"],
                                                                  "%sIndicesOther%s"%obj["muon"],
                                                                  ]]+[
            steps.jet.forwardFailedJetVeto( obj["jet"], ptAbove = 50, etaAbove = 3.5),
            steps.jet.uniquelyMatchedNonisoMuons(obj["jet"]),
            ])

    @staticmethod
    def selectionSteps(pars, withPlots = True) :
        obj = pars["objects"]
        bVar = ("%s"+pars["bVar"]+"%s")%calculables.jet.xcStrip(obj["jet"])
        lepton = obj[pars["lepton"]["name"]]
        lPtMin = pars["lepton"]["ptMin"]
        lEtaMax = pars["lepton"]["etaMax"]
        lIsoIndices = ("%s"+pars["selection"]["lIso"]["indices"]+"%s")%lepton

        topTag = pars['tag'].replace("Wlv","top").replace("QCD","top")
        selections = (
            [supy.steps.histos.multiplicity("%sIndices%s"%obj["jet"]),
             supy.steps.filters.multiplicity("%sIndices%s"%obj["jet"], **pars["nJets"]),
             
             supy.steps.histos.pt("mixedSumP4",100,0,300),
             supy.steps.filters.pt("mixedSumP4",min=20),
             
             topAsymmShell.lepIso(1,pars),
             supy.steps.filters.multiplicity("%sIndices%s"%lepton, max = 1), # drell-yann rejection
             topAsymmShell.lepIso(0,pars),

             supy.steps.filters.value(("%s"+pars["lepton"]["isoVar"]+"%s")%lepton, max = 1.0, indices = lIsoIndices, index = 0),
             supy.steps.filters.multiplicity("%sIndices%s"%lepton, min=pars["selection"]["lIso"]["N"],max=pars["selection"]["lIso"]["N"]),
             supy.steps.filters.pt("%sP4%s"%lepton, min = lPtMin, indices = lIsoIndices, index = 0),
             supy.steps.filters.absEta("%sP4%s"%lepton, max = lEtaMax, indices = lIsoIndices, index = 0),
             
             ]+[supy.steps.histos.value(bVar, 60,0,15, indices = "%sIndicesBtagged%s"%obj["jet"], index = i) for i in range(3)]+[
            calculables.jet.ProbabilityGivenBQN(obj["jet"], pars['bVar'], binning=(64,-1,15), samples = pars['topBsamples'], tag = topTag),
            supy.steps.histos.value("TopRatherThanWProbability", 100,0,1),
            #supy.steps.filters.value("TopRatherThanWProbability", min = 0.2),
            supy.steps.filters.value(bVar, indices = "%sIndicesBtagged%s"%obj["jet"], index = 1, min = 0.0),
            supy.steps.filters.value(bVar, indices = "%sIndicesBtagged%s"%obj["jet"], **pars["selection"]["bCut"]),
            ])
        return [s for s in selections if withPlots or s.isSelector or issubclass(type(s),supy.calculables.secondary)]

    @staticmethod
    def lepIso(index,pars) :
        lepton = pars["objects"][pars["lepton"]["name"]]
        return supy.steps.histos.value(("%s"+pars["lepton"]["isoVar"]+"%s")%lepton, 55,0,1.1, indices = "%sIndicesAnyIsoIsoOrder%s"%lepton, index=index)
