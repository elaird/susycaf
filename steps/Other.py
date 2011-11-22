import array,collections, ROOT as r
from supy import analysisStep,utils,steps
#####################################
class handleChecker(analysisStep) :
    def __init__(self, matches = ["handle", "Handle", "valid", "Valid"]) :
        self.run = False
        self.matches = matches

    def ofInterest(self, eventVars) :
        out = []
        for var in eventVars :
            for item in self.matches :
                if item in var :
                    out.append(var)
        return list(set(out))
    
    def uponAcceptance(self, eventVars) :
        if self.run : return
        self.run = True

        true = []
        false = []
        for var in self.ofInterest(eventVars) :
            if eventVars[var] :
                true.append(var)
            else :
                false.append(var)
        print "True:",sorted(true)
        print
        print "False:",sorted(false)
#####################################
class jsonMaker(analysisStep) :

    def __init__(self, calculateLumi = True) :
        self.lumisByRun = collections.defaultdict(list)
        self.calculateLumi = calculateLumi
        self.moreName="see below"

    def uponAcceptance(self,eventVars) :
        self.lumisByRun[eventVars["run"]].append(eventVars["lumiSection"])
    
    def varsToPickle(self) : return ["lumisByRun"]

    def outputSuffix(self) : return ".json"
    
    def mergeFunc(self, products) :
        d = collections.defaultdict(list)
        for lumisByRun in products["lumisByRun"] :
            for run,lumis in lumisByRun.iteritems() :
                d[run] += lumis

        d2 = {}
        for run,lumis in d.iteritems() :
            d2[run] = sorted(set(lumis))
            for ls in d2[run] :
                if 1 < lumis.count(ls) :
                    print "Run %d ls %d appears %d times in the lumiTree."%(run,ls,lumis.count(ls))

        json = utils.jsonFromRunDict(d2)
        lumi = utils.luminosity.recordedInvMicrobarns(json)/1e6 if self.calculateLumi else -1.0
        with open(self.outputFileName,"w") as file: print >> file, str(json).replace("'",'"')
        print "Wrote %.4f/pb json to : %s"%(lumi,self.outputFileName)
        print utils.hyphens
#####################################
class duplicateEventCheck(analysisStep) :
    def __init__(self) :
        self.events = collections.defaultdict(list)

    def uponAcceptance(self,ev) :
        self.events[(ev["run"], ev["lumiSection"])].append(ev["event"])

    def varsToPickle(self) : return ["events"]

    def mergeFunc(self, products) :
        def mergedEventDicts(l) :
            out = collections.defaultdict(list)
            for d in l :
                for key,value in d.iteritems() :
                    out[key] += value
            return out

        def duplicates(l) :
            s = set(l)
            for item in s :
                l.remove(item)
            return list(set(l))

        anyDups = False
        events = mergedEventDicts(products["events"])
        for runLs in sorted(events.keys()) :
            d = duplicates(events[runLs])
            if d :
                print "DUPLICATE EVENTS FOUND in run %d ls %d: %s"%(runLs[0], runLs[1], d)
                anyDups = True
        if not anyDups :
            print "No duplicate events were found."
#####################################






#####################################
class productGreaterFilter(analysisStep) :

    def __init__(self, threshold, variables, suffix = ""):
        self.threshold = threshold
        self.variables = variables
        self.moreName = "%s>=%.3f %s" % ("*".join(variables),threshold,suffix)

    def select (self,eventVars) :
        product = 1
        for var in self.variables : product *= eventVars[var]
        return product >= self.threshold
#####################################
class objectEtaSelector(analysisStep) :

    def __init__(self, cs, etaThreshold, index, p4String):
        self.index = index
        self.etaThreshold = etaThreshold
        self.cs = cs
        self.indicesName = "%sIndices%s" % self.cs
        self.p4sName = "%s%s%s" % (self.cs[0], p4String, self.cs[1])
        self.moreName = "%s%s; |eta[index[%d]]|<=%.1f" % (self.cs[0], self.cs[1], index, etaThreshold)

    def select (self,eventVars) :
        indices = eventVars[self.indicesName]
        if len(indices) <= self.index : return False
        p4s = eventVars[self.p4sName]
        return self.etaThreshold > abs(p4s.at(indices[self.index]).eta())
#####################################
class ptRatioLessThanSelector(analysisStep) :

    def __init__(self,numVar = None, denVar = None, threshold = None):
        for item in ["numVar","denVar","threshold"] :
            setattr(self,item,eval(item))
        self.moreName = "%s.pt() / %s.pt() < %.3f" % (numVar,denVar,threshold)

    def select (self,ev) :
        value = ev[self.numVar].pt() / ev[self.denVar].pt()
        return value<self.threshold
#####################################
class ptRatioHistogrammer(analysisStep) :

    def __init__(self,numVar = None, denVar = None):
        for item in ["numVar","denVar"] :
            setattr(self,item,eval(item))
        self.moreName = "%s.pt() / %s.pt()" % (self.numVar,self.denVar)

    def uponAcceptance (self,ev) :
        value = ev[self.numVar].pt() / ev[self.denVar].pt()
        self.book.fill(value,"ptRatio", 50, 0.0, 2.0, title = ";%s / %s;events / bin"%(self.numVar,self.denVar) )
#####################################
class runHistogrammer(analysisStep) :

    def __init__(self) :
        self.runDict={}
        
    def select (self,eventVars) :
        run = eventVars["run"]
        if run in self.runDict :
            self.runDict[run] += 1
        else :
            self.runDict[run] = 1
        return True

    def printStatistics(self) :
        printList=[]
        for key in sorted(self.runDict) :
            printList.append( (key,self.runDict[key]) )
        print self.name1(),printList
#####################################
class bxFilter(analysisStep) :

    def __init__(self,bxList) :
        self.bxList = bxList
        self.moreName = "[%s]" % ",".join(bxList)

    def select (self,eventVars) :
        return eventVars["bunch"] in self.bxList
#####################################
class pickEventSpecMaker(analysisStep) :
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/WorkBookPickEvents

    def __init__(self) :
        self.events = []
        
    def outputSuffix(self) :
        return "_pickEvents.txt"
    
    def uponAcceptance(self, eventVars) :
        self.events.append( (eventVars["run"], eventVars["lumiSection"], eventVars["event"]) )
        
    def varsToPickle(self) :
        return ["events"]

    def mergeFunc(self, products) :
        out = open(self.outputFileName, "w")
        for events in products["events"] :
            for event in events :
                out.write("%14d:%6d:%14d\n"%event)
        out.close()
        print "The pick events spec. file %s has been written."%self.outputFileName
class deadEcalFilter(analysisStep) :
    def __init__(self, jets = None, extraName = "", dR = None, dPhiStarCut = None) :
        for item in ["jets","extraName","dR","dPhiStarCut"] :
            setattr(self,item,eval(item))
        self.dps = "%sDeltaPhiStar%s%s"%(self.jets[0], self.jets[1], self.extraName)
        self.deadEcalDR = "%sDeadEcalDR%s%s"%(self.jets[0], self.jets[1], self.extraName)
        self.moreName = "%s%s; dR>%5.3f when deltaPhiStar<%5.3f"%(self.jets[0], self.jets[1], self.dR, self.dPhiStarCut)
        
    def select(self, eventVars) :
        passFilter = eventVars[self.dps][0][0]>self.dPhiStarCut or eventVars[self.deadEcalDR][0]>self.dR
        #if not passFilter :
        #    jet = eventVars["%sCorrectedP4%s"%self.jets].at(eventVars[self.dps]["DeltaPhiStarJetIndex"])
        #    self.book.fill(jet.pt(), "ptOfJetCausingDPhiVeto", 10, 0.0, 100.0, title = ";p_{T} of jet causing #Delta#phi* veto (GeV);events / bin")
        return passFilter
#####################################
class vertexHistogrammer(analysisStep) :
    def uponAcceptance(self,eV) :
        index = eV["vertexIndices"]
        sump3 = eV["vertexSumP3"]
        sumpt = eV["vertexSumPt"]
        if not len(index) : return
        
        #coord = reduce(lambda v,u: (v[0]+u[0],v[1]+u[1],v[2]+u[2]), [(sump3[i].x(),sump3[i].y(),sump3[i].z()) for i in index][1:], (0,0,0))
        #sump3Secondaries = type(sump3[0])(coord[0],coord[1],coord[2])

        self.book.fill( sumpt[index[0]], "vertex0SumPt", 40, 0, 1200, title = ";primary vertex #Sigma p_{T} (GeV); events / bin")
        self.book.fill( sump3[index[0]].rho(), "vertex0MPT%d", 40, 0, 400, title = ";primary vertex MPT (GeV);events / bin")

        self.book.fill( sum(map(sumpt.__getitem__,index[1:])), "vertexGt0SumPt", 100, 0, 400, title = ";secondary vertices #Sigma p_{T};events / bin")
        #self.book.fill( (sumpt[index[0]], sum(map(sumpt.__getitem__,index[1:]))), "vertexSumPt_0_all", (100,100), (0,0), (1200,400), title = ";primary vertex #Sigma p_{T};secondary vertices #Sigma p_{T};events / bin")
        #self.book.fill( (sump3[index[0]].rho(), sump3Secondaries.rho()), "vertexMPT_0_all", (100,100), (0,0), (400,200), title = ";primary vertex MPT;secondary verticies MPT;events / bin")

#####################################
class cutSorter(analysisStep) :
    def __init__(self, listOfSteps = [], applySelections = True ) :
        self.selectors = filter(lambda s: s.isSelector, listOfSteps)
        self.applySelections = applySelections
        self.moreName = "Applied" if applySelections else "Not Applied"
        self.bins = 1 << len(self.selectors)

    def select(self,eventVars) :
        selections = [s.select(eventVars) for s in self.selectors]
        self.book.fill( utils.intFromBits(selections), "cutSorterConfigurationCounts", 
                                  self.bins, -0.5, self.bins-0.5, title = ";cutConfiguration;events / bin")
        return (not self.applySelections) or all(selections)
        
    def endFunc(self, chains) :
        bins = len(self.selectors)
        self.book.fill(1, "cutSorterNames", bins, 0, bins, title = ";cutName", xAxisLabels = [sel.__class__.__name__ for sel in self.selectors])
        self.book.fill(1, "cutSorterMoreNames", bins, 0, bins, title = ";cutMoreName", xAxisLabels = [sel.moreName for sel in self.selectors])
#####################################
class compareMissing(analysisStep) :
    def __init__(self, missings, bins = 50, min = 0, max = 500) :
        self.missings = missings
        self.missingPairs = [(a,b) for a in missings for b in missings[missings.index(a)+1:]]
        for item in ["bins","min","max"] :
            setattr(self,item,eval(item))
            setattr(self,"%sPair"%item,(eval(item),eval(item)))
    def uponAcceptance(self,eV) :
        for m in self.missings :
            self.book.fill(eV[m].pt(), "%s.pt"%m, self.bins,self.min,self.max, title = ";%s.pt;events / bin"%m)
        for pair in self.missingPairs :
            self.book.fill((eV[pair[0]].pt(),eV[pair[1]].pt()), "%s.pt_vs_%s.pt"%pair, self.binsPair, self.minPair, self.maxPair,
                           title = ";%s.pt;%s.pt;events / bin"%pair)
#####################################
class smsMedianHistogrammer(analysisStep) :
    def __init__(self, cs) :
        self.cs = cs

        self.nBinsX = 12
        self.xLo =  400.0
        self.xHi = 1000.0

        self.nBinsY = 36
        self.yLo =  100.0
        self.yHi = 1000.0

        self.nBinsZ = 100
        self.zLo =    0.0
        self.zHi = 1000.0

    def nEvents(self, eventVars) :
        self.book.fill((eventVars["susyScanM0"], eventVars["susyScanM12"]), "nEvents",
                       (self.nBinsX, self.nBinsY), (self.xLo, self.yLo), (self.xHi, self.yHi),
                       title = ";m_{mother} (GeV);m_{LSP} (GeV);N events")

    def ht(self, eventVars) :
        var = "%sSumEt%s"%self.cs
        value = eventVars[var] if eventVars[var] else 0.0
        self.book.fill((eventVars["susyScanM0"], eventVars["susyScanM12"], value), var,
                       (self.nBinsX, self.nBinsY, self.nBinsZ), (self.xLo, self.yLo, self.zLo), (self.xHi, self.yHi, self.zHi),
                       title = ";m_{mother} (GeV);m_{LSP} (GeV);%s"%var)

    def mht(self, eventVars) :
        var = "%sSumP4%s"%self.cs
        value = eventVars[var].pt() if eventVars[var] else 0.0
        self.book.fill((eventVars["susyScanM0"], eventVars["susyScanM12"], value), var,
                       (self.nBinsX, self.nBinsY, self.nBinsZ), (self.xLo, self.yLo, self.zLo), (self.xHi, self.yHi, self.zHi),
                       title = ";m_{mother} (GeV);m_{LSP} (GeV);%s"%var)
        

    def jets(self, eventVars) :
        jets = eventVars["%sCorrectedP4%s"%self.cs] if eventVars["%sCorrectedP4%s"%self.cs] else []
        for i in range(2) :
            var = "%sJet%dPt%s"%(self.cs[0], i, self.cs[1])
            if len(jets)<i+1 : value = 0.0
            else :             value = jets.at(i).pt()
            self.book.fill((eventVars["susyScanM0"], eventVars["susyScanM12"], value), var,
                           (self.nBinsX, self.nBinsY, self.nBinsZ), (self.xLo, self.yLo, self.zLo), (self.xHi, self.yHi, self.zHi),
                           title = ";m_{mother} (GeV);m_{LSP} (GeV);%s"%var)

    def forwardJets(self, eventVars) :
        jets = eventVars["%sCorrectedP4%s"%self.cs] if eventVars["%sCorrectedP4%s"%self.cs] else []        
        forwardJets = filter(lambda x:abs(x.eta())>3.0, jets)
        value = max([jet.pt() for jet in forwardJets]) if forwardJets else 0.0
        var = "%sMaxForwardJetPt%s"%self.cs
        self.book.fill((eventVars["susyScanM0"], eventVars["susyScanM12"], value), var,
                       (self.nBinsX, self.nBinsY, self.nBinsZ), (self.xLo, self.yLo, self.zLo), (self.xHi, self.yHi, self.zHi),
                       title = ";m_{mother} (GeV);m_{LSP} (GeV);%s"%var)


    def uponAcceptance(self, eventVars) :
        self.nEvents(eventVars)
        self.ht(eventVars)
        self.mht(eventVars)
        self.jets(eventVars)
        self.forwardJets(eventVars)
        
    def outputSuffix(self) : return steps.master.outputSuffix()

    def oneHisto(self, name, zAxisTitle) :
        f = r.TFile(self.outputFileName, "UPDATE")
        h = f.Get("master/orFilter/%s"%name)
        outName = "%s_median"%name
        out = r.TH2D(outName, h.GetTitle(),
                     h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax(),
                     h.GetNbinsY(), h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax())
        out.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
        out.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
        out.GetZaxis().SetTitle(zAxisTitle)
        oneD = r.TH1D("%s_oneD"%name, h.GetTitle(), h.GetNbinsZ(), h.GetZaxis().GetXmin(), h.GetZaxis().GetXmax())
        
        for iBinX in range(1, 1+h.GetNbinsX()) :
            for iBinY in range(1, 1+h.GetNbinsY()) :
                oneD.Reset()
                for iBinZ in range(1, 1+h.GetNbinsZ()) :
                    oneD.SetBinContent(iBinZ, h.GetBinContent(iBinX, iBinY, iBinZ))
                if not oneD.Integral() : continue
                probSum = array.array('d', [0.5])
                q = array.array('d', [0.0]*len(probSum))
                oneD.GetQuantiles(len(probSum), q, probSum)
                out.SetBinContent(iBinX, iBinY, q[0])

        f.cd("master/orFilter")
        out.Write()
        r.gROOT.cd()
        f.Close()
        print "Output updated with %s."%name
        
    def mergeFunc(self, products) :
        self.oneHisto("%sSumEt%s"%self.cs, "Median HT (GeV)")
        self.oneHisto("%sSumP4%s"%self.cs, "Median MHT (GeV)")
        self.oneHisto("%sMaxForwardJetPt%s"%self.cs, "Median pT of (highest pT jet with |#eta| > 3.0) (GeV)")
        for i in range(2) :
            self.oneHisto("%sJet%dPt%s"%(self.cs[0], i, self.cs[1]), "Median pT of jet %d (GeV)"%(i+1))
#####################################


###   Obsolete   ####
#####################################
class variablePtGreaterFilter(analysisStep) :
    '''Obsolete: use ptFilter'''

    def __init__(self, threshold, variable, suffix = ""):
        self.threshold = threshold
        self.variable = variable
        self.moreName = "%s.pt()>=%.1f %s" % (variable,threshold,suffix)

    def select (self,eventVars) :
        return eventVars[self.variable].pt()>=self.threshold
#####################################
class variablePtLessFilter(analysisStep) :
    '''Obsolete: use ptFilter'''

    def __init__(self, threshold, variable, suffix = ""):
        self.threshold = threshold
        self.variable = variable
        self.moreName = "%s.pt()<=%.1f %s" % (variable,threshold,suffix)

    def select (self,eventVars) :
        return eventVars[self.variable].pt()<=self.threshold
#####################################
class vertexRequirementFilter(analysisStep) :
    '''Obsolete: use vertexIndices + multiplicityFilter'''

    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes#Good_Vertex_selection
    def __init__(self, minNdof = 5.0, maxAbsZ = 24.0, maxD0 = 2.0) :
        for item in ["minNdof","maxAbsZ","maxD0"]: setattr(self,item,eval(item))
        self.moreName = "any v: !fake; ndf>=%.1f; |z|<=%.1f cm; d0<=%.1f cm" % (minNdof,maxAbsZ,maxD0)

    def select(self,eventVars) :
        fake,ndof,pos = eventVars["vertexIsFake"], eventVars["vertexNdof"], eventVars["vertexPosition"]
        
        for i in range(pos.size()) :
            if fake.at(i) : continue
            if ndof.at(i) < self.minNdof : continue
            if abs(pos.at(i).Z()) > self.maxAbsZ : continue
            if abs(pos.at(i).Rho()) > self.maxD0 : continue
            return True
        return False
#####################################
class variableGreaterFilter(analysisStep) :
    '''Obsolete: use valueFilter'''

    def __init__(self, threshold, variable, suffix = ""):
        self.threshold = threshold
        self.variable = variable
        self.moreName = "%s>=%.3f %s" % (variable,threshold,suffix)

    def select (self,eventVars) :
        return eventVars[self.variable]>=self.threshold
#####################################
class variableLessFilter(analysisStep) :
    '''Obsolete: use valueFilter'''

    def __init__(self, threshold, variable, suffix = ""):
        self.threshold = threshold
        self.variable = variable
        self.moreName = "%s<%.3f %s" % (variable,threshold,suffix)

    def select (self,eventVars) :
        return eventVars[self.variable]<self.threshold
#####################################
