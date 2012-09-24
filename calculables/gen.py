from supy import wrappedChain,utils,calculables
import ROOT as r
##############################
class genSumP4(wrappedChain.calculable) :
    def update(self,_) :
        genP4 = self.source['genP4']
        self.value = genP4.at(4) + genP4.at(5)
##############################
class wNonQQbar(wrappedChain.calculable) :
    def update(self,_) :
        self.value = None if self.source['genQQbar'] else 1
##############################
class wQQbar(wrappedChain.calculable) :
    def update(self,_) :
        self.value = 1 if self.source['genQQbar'] else None
##############################
class genQQbar(wrappedChain.calculable) :
    def update(self,_) :
        if self.source['isRealData'] :
            self.value = ()
            return
        ids = list(self.source['genPdgId'])
        iHard = self.source['genIndicesHardPartons']
        self.value = tuple(sorted(iHard,key = ids.__getitem__,reverse = True)) \
                     if not sum([ids[i] for i in iHard]) else tuple()
##############################
class genIndicesHardPartons(wrappedChain.calculable) :
    def __init__(self,indices = (4,5)) : self.value = indices
    def update(self,_) : pass
##############################
class genStatus1P4(wrappedChain.calculable) :
    def update(self,_) :
        self.value = []
        for i in range(self.source["genP4"].size()) :
            if self.source["genStatus"].at(i)!=1 : continue
            self.value.append(self.source["genP4"][i])
##############################
class genIndices(wrappedChain.calculable) :
    @property
    def name(self) : return "genIndices" + self.label

    def __init__(self, pdgs = [], label = None, status = [], motherPdgs = []) :
        self.label = label
        self.PDGs = frozenset(pdgs)
        self.status = frozenset(status)
        self.motherPdgs = frozenset(motherPdgs)
        self.moreName = "; ".join(["pdgId in %s" %str(list(self.PDGs)),
                                   "status in %s"%str(list(self.status)),
                                   "motherPdg in %s"%str(list(self.motherPdgs))
                                   ])

    def update(self,_) :
        pdg = self.source["genPdgId"]
        status = self.source["genStatus"]
        motherPdg = self.source["genMotherPdgId"]
        self.value = filter( lambda i: ( (not self.PDGs) or (pdg.at(i) in self.PDGs) ) and \
                                 ( (not self.status) or (status.at(i) in self.status) ) and \
                                 ( (not self.motherPdgs) or (motherPdg.at(i) in self.motherPdgs) ),
                             range(pdg.size()) )

class genIndicesPtSorted(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sPtSorted"%self.label

    def __init__(self, label = "") :
        self.label = "genIndices"+label

    def update(self,_) :
        p4 = self.source["genP4"]
        self.value = sorted(self.source[self.label], key = lambda i:p4.at(i).pt(), reverse = True)

class genRootSHat(wrappedChain.calculable) :
    def update(self,_) :
        iHard = self.source["genIndicesHardPartons"]
        p4s = self.source["genP4"]
        self.value = None if not iHard else (p4s.at(iHard[0])+p4s.at(iHard[1])).mass()

class genSumPt(wrappedChain.calculable) :
    @property
    def name(self) :
        return "_".join(["genSumPt"]+self.indexLabels)

    def __init__(self, indexLabels = []) :
        self.indexLabels = map(lambda s:s.replace("genIndices",""), indexLabels)

    def update(self,_) :
        indices = []
        for label in self.indexLabels :
            indices += self.source["genIndices"+label]
        indices = set(indices)

        self.value = 0.0
        p4 = self.source["genP4"]
        for i in indices :
            self.value += p4.at(i).pt()

class DeltaR(wrappedChain.calculable) :
    def __init__(self, indices = "") :
        self.fixes = (indices,"")

    def update(self,_) :
        p4 = self.source["genP4"]
        indices = self.source[self.fixes[0]]
        self.value = None
        if len(indices)==2 :
            self.value = r.Math.VectorUtil.DeltaR(p4.at(indices[0]), p4.at(indices[1]))

class MinDeltaPhiMet(wrappedChain.calculable) :
    def __init__(self, indices = "", met = "") :
        self.fixes = (indices, met)

    def update(self,_) :
        p4 = self.source["genP4"]
        indices = self.source[self.fixes[0]]
        met = self.source[self.fixes[1]]
        self.value = min([abs(r.Math.VectorUtil.DeltaPhi(met, p4.at(i))) for i in indices])

class JetIndices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, etaMax = None) :
        self.fixes = collection
        self.p4 = "gen%sJetsP4"%collection[0]
        self.ptMin = ptMin
        self.etaMax = etaMax
        self.moreName = "pT>%g GeV; |eta|<%g"%(self.ptMin, self.etaMax)

    def update(self,_) :
        p4 = self.source[self.p4]
        self.value = []
        for i in range(p4.size()) :
            jet = p4.at(i)
            if jet.pt()<self.ptMin : continue
            if abs(jet.eta())>self.etaMax : continue
            self.value.append(i)
##############################
class susyIniIndices(wrappedChain.calculable) :
    def __init__(self) :
        self.moreName = "status 3; SUSY; non-SUSY parents"

    def update(self,_) :
        self.value = []
        if not self.source["genHandleValid"] : return
        nParticles = len(self.source["genPdgId"])
        for iParticle in range(nParticles) :
            if self.source["genStatus"].at(iParticle)!=3 : #consider only status 3 particles
                continue
            if not utils.isSusy(self.source["genPdgId"].at(iParticle)) : #which are SUSY particles
                continue
            if utils.isSusy(self.source["genMotherPdgId"].at(iParticle)) : #whose mothers are not SUSY particles
                continue
            self.value.append(iParticle)
##############################
class genIndicesStatus3Status3SusyMother(wrappedChain.calculable) :
    def __init__(self) :
        self.moreName = "status 3; non-SUSY; SUSY parents"

    def update(self,_) :
        self.value = []
        if not self.source["genHandleValid"] : return
        nParticles = len(self.source["genPdgId"])
        for iParticle in range(nParticles) :
            if self.source["genStatus"].at(iParticle)!=3 : #consider only status 3 particles
                continue
            if utils.isSusy(self.source["genPdgId"].at(iParticle)) : #which are not SUSY particles
                continue
            if not utils.isSusy(self.source["genMotherPdgId"].at(iParticle)) : #whose mothers are SUSY particles
                continue
            self.value.append(iParticle)
##############################
#class nonSusyFromSusySumEt(wrappedChain.calculable) :
#    def __init__(self) : #, collection = None) :
##        self.fixes = collection
#        self.stash(["genIndicesStatus3Status3SusyMother","P4"])
#    def update(self, ignored) :
#        p4s = self.source[self.P4]
#        indices = self.source["genIndicesStatus3Status3SusyMother"]
#        self.value = reduce( lambda x,i: x+p4s.at(i).Et(), indices , 0)
###############################
class nonSusyFromSusySumP4(wrappedChain.calculable) :
    def update(self,_) :
        indices = self.source["genIndicesStatus3Status3SusyMother"]
        assert len(indices)==4,indices
        self.value = utils.LorentzV()
        for i in indices :
            self.value += self.source["genP4"].at(i)
##############################
class susyIniSumP4(wrappedChain.calculable) :
    def update(self,_) :
        indices = self.source["susyIniIndices"]
        assert len(indices)==2,indices
        self.value = utils.LorentzV()
        for i in indices :
            self.value += self.source["genP4"].at(i)
##############################
class SumP4(wrappedChain.calculable) :
    def __init__(self, indices = "") :
        self.fixes = ("", indices)

    def update(self,_) :
        self.value = utils.LorentzV()
        for i in self.source[self.fixes[1]] :
            self.value += self.source["genP4"].at(i)
##############################
class genIndicesB(wrappedChain.calculable) :
    def update(self,_) :
        ids = self.source['genPdgId']
        self.value = filter(lambda i: abs(ids[i]) is 5, range(len(ids)))
##############################
class genIndicesWqq(wrappedChain.calculable) :
    def update(self,_) :
        ids = self.source['genPdgId']
        mom = self.source['genMotherPdgId']
        self.value = filter(lambda i: abs(mom[i]) is 24 and abs(ids[i]) < 5, range(len(ids)))
##############################
class genIndicesStatus3NoStatus3Daughter(wrappedChain.calculable) :
    def update(self,_) :
        status = self.source["genStatus"]
        mother = self.source["genMotherIndex"]

        status3List = filter( lambda i: status.at(i)==3, range(status.size()) )
        motherIndices = set([mother[i] for i in status3List])
        self.value = filter( lambda i: i not in motherIndices, status3List )
##############################
class genMinDeltaRPhotonOther(wrappedChain.calculable) :
    @property
    def name(self) : return "genMinDeltaRPhotonOther"+self.label
    
    def __init__(self, label) :
        self.label = label
        
    def update(self,_) :
        st3Indices = self.source["genIndicesStatus3NoStatus3Daughter"]
        genP4s = self.source["genP4"]
        ids = self.source["genPdgId"]

        def minDeltaR(photonIndex) :
            candidates = [ (r.Math.VectorUtil.DeltaR(genP4s.at(photonIndex), genP4s.at(i)), ids.at(i))  for i in st3Indices]
            return min(filter(lambda x: x[1]!=22, candidates))

        photonIndices = self.source["genIndices"+self.label]
        self.value = min( [minDeltaR(photonIndex) for photonIndex in photonIndices] )[0] if len(photonIndices) else None
##############################
class genIsolations(wrappedChain.calculable) :
    @property
    def name(self) : return "genIsolation"+self.label
    
    def __init__(self, label = None, coneSize = None) :
        for item in ["label","coneSize"] :
            setattr(self,item,eval(item))
        
    def update(self, _) :
        genP4s = self.source["genP4"]
        nGen = genP4s.size()
        genIndices = self.source["genIndices"+self.label]
        self.value = {}
        for genIndex in genIndices :
            iso = 0.0
            for iParticle in range(nGen) :
                if iParticle==genIndex : continue
                if self.source["genStatus"].at(iParticle)!=1 : continue
                if self.source["genPdgId"].at(iParticle) in [-12, 12, -14, 14, -16, 16] : continue
                if r.Math.VectorUtil.DeltaR(genP4s.at(genIndex), genP4s.at(iParticle)) > self.coneSize : continue
                iso += genP4s.at(iParticle).pt()
            self.value[genIndex] = iso
##############################
class genPhotonCategory(wrappedChain.calculable) :
    @property
    def name(self) :
        return "category"+self.label

    def __init__(self, label) :
        self.label = label
        
    def update(self, _) :
        self.value = {}

        for index in self.source["genIndices"+self.label] :
            moId = self.source["genMotherPdgId"][index]
            if moId==22 :
                self.value[index] = "photonMother"
            elif abs(moId)<22 :
                self.value[index] = "quarkMother"
            else :
                self.value[index] = "otherMother"
##############################
class genParticleCounter(wrappedChain.calculable) :
    @property
    def name(self) : return "GenParticleCategoryCounts"

    def __init__(self) :
        self.value = {}
        self.pdgToCategory = {}

        #copied from PDG
        self.initPdgToCategory( 1, 6,"quark")
        self.initPdgToCategory(21,21,"gluon")

        self.initPdgToCategory(1000001,1000004,"squarkL")#left-handed
        self.initPdgToCategory(1000005,1000006,"squarkA")#ambiguous
        self.initPdgToCategory(1000011,1000016,"slepton")
        self.initPdgToCategory(2000001,2000004,"squarkR")#right-handed
        self.initPdgToCategory(2000005,2000006,"squarkA")#ambiguous
        self.initPdgToCategory(2000011,2000011,"slepton")
        self.initPdgToCategory(2000013,2000013,"slepton")
        self.initPdgToCategory(2000015,2000015,"slepton")
        self.initPdgToCategory(1000021,1000021,"gluino")
        self.initPdgToCategory(1000022,1000023,"chi0")
        self.initPdgToCategory(1000024,1000024,"chi+")
        self.initPdgToCategory(1000025,1000025,"chi0")
        self.initPdgToCategory(1000035,1000035,"chi0")
        self.initPdgToCategory(1000037,1000037,"chi+")
        self.initPdgToCategory(1000039,1000039,"gravitino")

        self.combineCategories(["squarkL","squarkR","squarkA"], "squark")
        self.combineCategories(["slepton","chi0","chi+","gravitino"], "otherSusy")

        self.badCategoryName = "noName"
        self.categories = list(set(self.pdgToCategory.values()))
        self.categories.append(self.badCategoryName)
        self.categories.sort()
        #self.printDict(self.pdgToCategory)

    def initPdgToCategory(self,lower,upper,label) :
        for i in range(lower,upper+1) :
            self.pdgToCategory[i]=label
        for i in range(-upper,-lower+1) :
            self.pdgToCategory[i]=label

    def combineCategories(self,someList,someLabel) :
        for key in self.pdgToCategory :
            if self.pdgToCategory[key] in someList :
                self.pdgToCategory[key]=someLabel
        
    def printDict(self,someDict) :
        for key in someDict :
            print key,someDict[key]

    def zeroCategoryCounts(self) :
        for key in self.categories :
            self.value[key]=0

    def incrementCategory(self,pdgId) :
        if pdgId in self.pdgToCategory:
            category=self.pdgToCategory[pdgId]
        else :
            category=self.badCategoryName
        self.value[category]+=1
        #print "found one:",iParticle,pdgId

    def update(self,_) :
        self.zeroCategoryCounts()
        if not self.source["genHandleValid"] : return
        nParticles = len(self.source["genPdgId"])

        #Susy counts
        for iParticle in self.source["susyIniIndices"] :
            self.incrementCategory(self.source["genPdgId"].at(iParticle))

        #initial state counts
        for iParticle in range(nParticles) :
            #consider only status 3 particles
            if self.source["genStatus"].at(iParticle)!=3 : continue
            #whose mothers are protons
            if self.source["genMotherPdgId"].at(iParticle)!=2212 : continue
            #whose mothers have index 0 or 1
            if self.source["genMotherIndex"].at(iParticle) not in [0,1] : continue
            self.incrementCategory(self.source["genPdgId"].at(iParticle))
######################################
class qDirExpectation(calculables.secondary) :
    var = ""
    limit = 1
    tag = ""
    sample = ""
    p = None

    def onlySamples(self) : return [self.sample]

    def setup(self,*_) :
        import numpy as np
        orig = self.fromCache( [self.sample], [self.var], tag = self.tag)[self.sample][self.var]
        if not orig :
            print "No cache: %s; %s"%(self.sample,str(self))
            return
        values = [orig.GetBinContent(i) for i in range(orig.GetNbinsX()+2)]
        neighbors = 10
        for i in range(neighbors,len(values)-neighbors) : orig.SetBinContent(i, sum(values[i-neighbors:i+neighbors+1]) / (1+2*neighbors))

        edges = utils.edgesRebinned(orig, targetUncRel = 0.1)
        hist = orig.Rebin(len(edges)-1, orig.GetName()+"_rebinned", edges)
        vals  = [hist.GetBinContent(i) for i in range(1,len(edges))]
        del hist
        iZero = edges.index(0)
        R = np.array(vals[iZero:])
        L = np.array(vals[:iZero])[::-1]
        p = (R-L) / ( R + L )

        self.p = r.TH1D(self.name, ";|%s|;|<qDir>|"%self.var, len(edges[iZero:])-1, edges[iZero:])
        for i in range(len(p)) : self.p.SetBinContent(i+1,p[i])
        self.p.SetBinContent(len(edges[iZero:])+2, edges[-1])

        widths = [high-low for low,high in zip(edges[iZero:-1],edges[iZero+1:])]
        q = (R+L)/(widths * (sum(R)+sum(L)))
        self.q = self.p.Clone(self.name+"_pdist")
        self.q.Reset()
        self.q.SetTitle(";|%s|;probability of |%s|"%(self.var,self.var))
        for i in range(len(q)) : self.q.SetBinContent(i+1,q[i])

        self.mean = sum(a*b for a,b in zip(p,(R+L))) / (sum(R)+sum(L))
        self.value = self.calculate

    def reportCache(self) :
        fileName = '/'.join(self.outputFileName.split('/')[:-1]+[self.name])
        self.setup()
        if not self.p : return
        c = r.TCanvas()
        self.p.SetLineWidth(2)
        self.p.SetMaximum(1)
        self.p.SetMinimum(0)
        self.p.Draw('hist')
        self.q.SetLineWidth(2)
        self.q.SetLineColor(r.kRed)
        self.q.Draw('hist same')
        mean = self.q.Clone('mean'); mean.Reset();
        for i in range(mean.GetNbinsX()) : mean.SetBinContent(i+1,self.mean)
        mean.SetTitle("%.2f"%self.mean)
        mean.SetLineColor(r.kBlue)
        mean.Draw('hist same')
        utils.tCanvasPrintPdf(c,fileName)
        del mean
        del c

    def update(self,_) : pass

    def calculate(self, top, tbar) :
        var = self.varFunction(top,tbar)
        p = max(0,self.p.GetBinContent(self.p.FindFixBin(abs(var)))) if self.p else 0
        return p if var > 0 else -p

    def uponAcceptance(self,ev) :
        if ev['isRealData'] : return
        qqbar = ev['genQQbar']
        if not qqbar : return
        iTT = ev['genTopTTbar']
        if not iTT : return
        p4 = ev['genP4']
        qdir = 1 if p4[qqbar[0]].pz()>0 else -1
        var = self.varFunction(p4[iTT[0]],p4[iTT[1]])
        self.book.fill(qdir*var, self.var, 1000, -self.limit, self.limit, title = ";qdir * %s;events / bin"%self.var )

    def varFunction(self,top,tbar) : return


class qDirExpectation_(qDirExpectation) :
    def __init__(self, var, limit, tag, sample) :
        for item in ['var','tag','sample', 'limit'] : setattr(self,item,eval(item))
        self.fixes = ('',var)

    def varFunction(self,top,tbar) : return self.source[self.var]

class qDirExpectation_SumRapidities(qDirExpectation) :
    def varFunction(self,top,tbar) : return top.Rapidity() + tbar.Rapidity()
    def __init__(self, tag, sample) :
        self.var = "SumRapidities"
        self.limit = 5
        self.tag = tag
        self.sample = sample

class qDirExpectation_EtaSum(qDirExpectation) :
    def varFunction(self,top,tbar) : return (top + tbar).Eta()
    def __init__(self, tag, sample) :
        self.var = "EtaSum"
        self.limit = 8
        self.tag = tag
        self.sample = sample

class qDirExpectation_RapiditySum(qDirExpectation) :
    def varFunction(self,top,tbar) : return (top + tbar).Rapidity()
    def __init__(self, tag, sample) :
        self.var = "RapiditySum"
        self.limit = 2.5
        self.tag = tag
        self.sample = sample

##############################
class isrWeight(wrappedChain.calculable) :
    def __init__(self, model = "", var = "") :
        assert model in ["T1", "T2"],model
        self.model = model
        self.var = var
        self.histos = {}

    def setup(self) :
        f = r.TFile("data/ISRWeights_Topology%s.root"%self.model)
        for item in f.GetListOfKeys() :
            name = item.GetName()
            h = f.Get(name).Clone()
            h.SetDirectory(0)
            assert name.startswith("h_ISRWeight_lastPt_")
            mX,mY = name.split("_")[-2:]
            self.histos[(float(mX),float(mY))] = h
        f.Close()

    def update(self,_) :
        if not self.histos :
            self.setup()

        h = self.histos[(self.source["susyScanmGL"], self.source["susyScanmLSP"])]
        x = self.source["susyIniSumP4"].pt()
        self.value = h.GetBinContent(h.FindBin(x))
