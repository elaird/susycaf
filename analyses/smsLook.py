#!/usr/bin/env python

#import os,analysis,steps,calculables,samples,organizer,plotter,utils
import os,supy,steps,calculables,samples
import ROOT as r

class smsLook(supy.analysis) :
    def parameters(self) :
        objects = self.vary()
        fields =                                                [ "jet",            "met",            "muon",        "electron",        "photon",
                                                                  "compJet",    "compMet",
                                                                  "rechit", "muonsInJets", "jetPtMin", "jetId"]

        objects["caloAK5JetMet_recoLepPhot"] = dict(zip(fields, [("xcak5Jet","Pat"),"metP4TypeIPF",("muon","Pat"),("electron","Pat"),("photon","Pat"),
                                                                 ("xcak5JetPF","Pat"),"metP4PF",
                                                                 "Calo",     False,         50.0,      "JetIDloose"]))
        
        return { "objects": objects,
                 "nJetsMinMax" :      self.vary(dict([ ("ge2",(2,None)),  ("2",(2,2)),  ("ge3",(3,None)) ]       [0:1] )),
                 "mcSoup" :           self.vary(dict([ ("pythia6","py6"), ("pythia8","py8"), ("madgraph","mg") ] [0:1] )),
                 "etRatherThanPt" : [True,False]        [0],
                 "lowPtThreshold" : 30.0,
                 "lowPtName" : "lowPt",
                 "triggerList" : ("HLT_HT100U","HLT_HT100U_v3","HLT_HT120U","HLT_HT140U","HLT_HT150U_v3"),#required to be a sorted tuple
                 }

    def calcListJet(self, obj, etRatherThanPt, lowPtThreshold, lowPtName) :
        def calcList(jet, met, photon, muon, electron, muonsInJets, jetPtMin, jetIdFlag) :
            outList = [
                calculables.xclean.xcJet(jet,
                                         applyResidualCorrectionsToData = False,
                                         gamma = photon,
                                         gammaDR = 0.5,
                                         muon = muon,
                                         muonDR = 0.5,
                                         correctForMuons = not muonsInJets,
                                         electron = electron,
                                         electronDR = 0.5),
                calculables.jet.Indices( jet, jetPtMin, etaMax = 3.0, flagName = jetIdFlag),
                calculables.jet.Indices( jet, lowPtThreshold, etaMax = 3.0, flagName = jetIdFlag, extraName = lowPtName),
                
                calculables.jet.SumP4(jet),
                calculables.jet.SumP4(jet, extraName = lowPtName),
                calculables.jet.DeltaPhiStar(jet, extraName = lowPtName),
                calculables.jet.DeltaPseudoJet(jet, etRatherThanPt),
                calculables.jet.AlphaT(jet, etRatherThanPt),
                calculables.jet.AlphaTMet(jet, etRatherThanPt, met),
                calculables.jet.MhtOverMet(jet, met),
                calculables.jet.deadEcalDR(jet, extraName = lowPtName, minNXtals = 10),
                ]
            return outList + supy.calculables.fromCollections(calculables.jet, [jet])

        outList = calcList(obj["jet"], obj["met"], obj["photon"], obj["muon"], obj["electron"], obj["muonsInJets"], obj["jetPtMin"], obj["jetId"])
        if obj["compJet"]!=None and obj["compMet"]!=None :
            outList += calcList(obj["compJet"], obj["compMet"], obj["photon"], obj["muon"], obj["electron"], obj["muonsInJets"], obj["jetPtMin"], obj["jetId"])
        return outList

    def calcListOther(self, obj, triggers) :
        return [
            calculables.xclean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),

            calculables.muon.Indices( obj["muon"], ptMin = 10, combinedRelIsoMax = 0.15),
            calculables.electron.Indices( obj["electron"], ptMin = 10, simpleEleID = "95", useCombinedIso = True),
            calculables.photon.Indices(obj["photon"], ptMin = 25, flagName = "photonIDLooseFromTwikiPat"),
            
            calculables.vertex.ID(),
            calculables.vertex.Indices(),
            #calculables.Other.lowestUnPrescaledTrigger(triggers),
            ]
    
    def listOfCalculables(self, params) :
        obj = params["objects"]
        outList  = supy.calculables.zeroArgs(supy.calculables)
        outList += supy.calculables.zeroArgs(calculables)
        outList += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])
        outList += supy.calculables.fromCollections(calculables.electron, [obj["electron"]])
        outList += supy.calculables.fromCollections(calculables.photon, [obj["photon"]])
        #outList += self.calcListOther(obj, params["triggerList"])
        outList += self.calcListJet(obj, params["etRatherThanPt"], params["lowPtThreshold"], params["lowPtName"])
        return outList
    
    def listOfSteps(self, params) :
        _jet  = params["objects"]["jet"]
        _electron = params["objects"]["electron"]
        _muon = params["objects"]["muon"]
        _photon = params["objects"]["photon"]
        _met  = params["objects"]["met"]
        _etRatherThanPt = params["etRatherThanPt"]
        _et = "Et" if _etRatherThanPt else "Pt"

        return [
            supy.steps.printer.progressPrinter(),
            supy.steps.other.collector(["SimpModelScanmGL","SimpModelScanmLSP"]),
            #supy.steps.other.orFilter([steps.Gen.ParticleCountFilter({"gluino":2}), steps.Gen.ParticleCountFilter({"squark":2})]),
            #steps.other.smsMedianHistogrammer(_jet),
            ]+[
            supy.steps.histos.generic(("SimpModelScanmGL","SimpModelScanmLSP"),(101,38),(i,i),(2020+i,760+i), title="%d;m_{0};m_{1/2};events/bin"%i) for i in range(20)]
    
    def listOfSampleDictionaries(self) :
        sampleDict = supy.samples.SampleHolder()
   
        sampleDict.add("t1_1000_50", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_50/t1_1000_50.root"]', lumi = 1.1e3)
        sampleDict.add("t1_1000_600", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim1000_600/t1_1000_600.root"]', lumi = 1.1e3)
        sampleDict.add("t1_400_300", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim400_300/t1_400_300.root"]', lumi = 1.1e3)        
        sampleDict.add("t1_3_points", '["/uscms/home/yeshaq/nobackup/supy-output/smsSkim/sms_3_points.root"]', lumi = 1.1e3)

        return [sampleDict]
        #return [samples.mc]

    def listOfSamples(self,params) :
        from supy.samples import specify
        #return specify( #nFilesMax = 10, nEventsMax = 20000,
        #               names = ["t1.yos"])
        return supy.samples.specify(names = ["t1_3_points"])
                                    #"t1_1000_50",
                                    #"t1_1000_600",
                                    #"t1_400_300"
                                    
   
    def conclude(self,pars) :
        org = self.organizer(pars)

        #plot
        pl = supy.plotter(org,
                             pdfFileName = self.pdfFileName(org.tag),
                             #samplesForRatios = ("2010 Data","standard_model"),
                             #sampleLabelsForRatios = ("data","s.m."),
                             #whiteList = ["lowestUnPrescaledTrigger"],
                             doLog = False,
                             #compactOutput = True,
                             #noSci = True,
                             #pegMinimum = 0.1,
                             blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                             )
        pl.plotAll()
