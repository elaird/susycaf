from supy import wrappedChain,utils,configuration
import math
try:
    import numpy as np
except:
    pass

class pt(wrappedChain.calculable) :
    def __init__(self,var) :
        self.fixes = ("%s."%var, "")
        self.var = var
    def update(self,_) : self.value = self.source[self.var].pt()
#####################################
class eta(wrappedChain.calculable) :
    def __init__(self,var) :
        self.fixes = ("%s."%var, "")
        self.var = var
    def update(self,_) : self.value = self.source[self.var].eta()
#####################################
class ecalDeadTowerIsBarrel(wrappedChain.calculable) :
    etaBE = configuration.detectorSpecs()["cms"]["etaBE"]
    def update(self,ignored) : self.value = map( self.isBarrel, self.source["ecalDeadTowerTrigPrimP4"] )
    def isBarrel(self, p4) : return abs(p4.eta()) < self.etaBE
#####################################
class lowestUnPrescaledTrigger(wrappedChain.calculable) :
    def __init__(self, sortedListOfPaths = []) :
        self.sortedListOfPaths = sortedListOfPaths
        self.cached = dict()
        self.moreName = "lowest unprescaled of "+','.join(self.sortedListOfPaths).replace("HLT_","")
        
    def update(self, ignored) :
        key = (self.source["run"],self.source["lumiSection"])
        if key not in self.cached :
            self.cached[key] = None
            for path in self.sortedListOfPaths :
                if self.source["prescaled"][path]==1 :
                    self.cached[key] = path
                    break
        self.value = self.cached[key]
##############################
class pthatLess(wrappedChain.calculable) :
    def __init__(self, maxPtHat = None) :
        self.fixes = ("","%d"%maxPtHat)
        self.maxPtHat = maxPtHat
    def update(self,ignored) : self.value = None if self.maxPtHat < self.source["genpthat"] else 1
##############################
class Mt(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sMt%s%s"%(self.fixes[0], self.fixes[1], self.met)
    
    def __init__(self, collection = None, met = None, byHand = True , allowNonIso = False, isSumP4=False) :
        self.met = met
        self.fixes = collection
        self.stash(["Indices","IndicesNonIso","P4"])
        self.byHand = byHand
        self.isSumP4 = isSumP4
        self.allowNonIso = allowNonIso
        self.moreName = "%s%s, %s, byHand=%d"%(collection[0], collection[1], met, byHand)

    def update(self, ignored) :
        if (not self.source[self.Indices]) and not (self.allowNonIso or self.source[self.indicesNonIso]) :
            self.value= -1.0
            return
        index = self.source[self.Indices][0] if self.source[self.Indices] else \
                self.source[self.IndicesNonIso][0]
        lep = self.source[self.P4][index]
        met = self.source[self.met] * (-1 if self.isSumP4 else 1)

        if self.byHand :
            self.value = math.sqrt( 2.0*lep.pt()*met.pt()*(1.0 - math.cos(r.Math.VectorUtil.DeltaPhi(lep, met))) )
        else :
            self.value = (lep+met).Mt()
##############################
class SumP4(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","P4"])
    def update(self, ignored) :
        p4s = self.source[self.P4]
        indices = self.source[self.Indices]
        self.value = reduce( lambda x,i: x+p4s.at(i), indices, utils.LorentzV()) if len(indices) else None
##############################
class SumEt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","P4"])
    def update(self, ignored) :
        p4s = self.source[self.P4]
        indices = self.source[self.Indices]
        self.value = reduce( lambda x,i: x+p4s.at(i).Et(), indices , 0)
##############################
class RecHitSumPt(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sRecHitSumPt"%self.collection
    def __init__(self, collection = None, minEcalSeverity = 2, minHcalSeverity = 10, excludeHf = True, prune = True) :
        for item in ["collection", "minEcalSeverity", "minHcalSeverity", "prune"] :
            setattr(self, item, eval(item))
        self.considerSeverity = self.collection=="Calo"
        self.subdetectors = configuration.detectorSpecs()["cms"]["%sSubdetectors"%self.collection]
        if excludeHf : self.subdetectors = filter(lambda s:"Hf" not in s, self.subdetectors)
        self.recHitCollections = configuration.detectorSpecs()["cms"]["%sRecHitCollections"%self.collection]
    def update(self, ignored) :
        self.value = 0.0
        for detector in self.subdetectors :
            minSeverityLevel = self.minEcalSeverity if detector[0]=="E" else self.minHcalSeverity
            considered = []
            for collectionName in self.recHitCollections :
                p4Var = "rechit%s%s%s%s"%(collectionName, self.collection, "P4", detector)
                slVar = "rechit%s%s%s%s"%(collectionName, self.collection, "SeverityLevel", detector)
                for iHit in range(len(self.source[p4Var])) :
                    if self.considerSeverity and (self.source[slVar].at(iHit)<minSeverityLevel) : continue
                    hit = self.source[p4Var].at(iHit)
                    coords = (hit.eta(), hit.phi())
                    if self.prune :
                        if coords in considered : continue
                        considered.append(coords)
                    self.value += hit.pt()
##############################
class RecHitSumP4(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sRecHitSumP4"%self.collection
    def __init__(self, collection = None, minEcalSeverity = 2, minHcalSeverity = 10) :
        for item in ["collection", "minEcalSeverity", "minHcalSeverity"] :
            setattr(self, item, eval(item))
        self.considerSeverity = self.collection=="Calo"
        self.subdetectors = configuration.detectorSpecs()["cms"]["%sSubdetectors"%self.collection]
        self.recHitCollections = configuration.detectorSpecs()["cms"]["%sRecHitCollections"%self.collection]
    def update(self, ignored) :
        self.value = utils.LorentzV()
        for detector in self.subdetectors :
            minSeverityLevel = self.minEcalSeverity if detector[0]=="E" else self.minHcalSeverity
            for collectionName in self.recHitCollections :
                p4Var = "rechit%s%s%s%s"%(collectionName, self.collection, "P4", detector)
                slVar = "rechit%s%s%s%s"%(collectionName, self.collection, "SeverityLevel", detector)
                for iHit in range(len(self.source[p4Var])) :
                    if self.considerSeverity and (self.source[slVar].at(iHit)<minSeverityLevel) : continue
                    self.value += self.source[p4Var].at(iHit)
##############################
class metPlusParticles(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sPlus%s%s"%(self.met, self.particles[0], self.particles[1])
    def __init__(self, met, particles) :
        self.met = met
        self.particles = particles
        self.moreName = "%s + %s%s"%(self.met, self.particles[0], self.particles[1])
    def update(self, ignored) :
        self.value = self.source[self.met] + self.source["%sSumP4%s"%self.particles]
##############################
class minDeltaRToJet(wrappedChain.calculable) :
    @property
    def name(self) : return "%s%sMinDeltaRToJet%s%s"% (self.particles[0], self.particles[1], self.jets[0], self.jets[1])

    def __init__(self, particles, jets) :
        for item in ["particles","jets"] :
            setattr(self,item,eval(item))
        self.particleIndices = "%sIndices%s"%self.particles
        self.particleP4s     = "%sP4%s"     %self.particles

        self.jetIndices = "%sIndices%s"    %self.jets
        self.jetP4s     = "%sCorrectedP4%s"%self.jets

    def update(self, ignored) :
        self.value = {}
        particleIndices = self.source[self.particleIndices]
        particles       = self.source[self.particleP4s]

        jetIndices    = self.source[self.jetIndices]
        jets          = self.source[self.jetP4s]
        for iParticle in particleIndices :
            self.value[iParticle] = min([r.Math.VectorUtil.DeltaR( particles.at(iParticle), jets.at(iJet) ) for iJet in jetIndices]) if len(jetIndices) else None
#####################################
class jsonWeight(wrappedChain.calculable) :
    def __init__(self, fileName = "", acceptFutureRuns = False) :
        self.moreName = "run:ls in %s"%fileName
        self.acceptFutureRuns = acceptFutureRuns
        if self.acceptFutureRuns : self.moreName += " OR future runs"

        self.json = {}
        self.runs = []
        self.maxRunInJson = -1
        if fileName :
            file = open(fileName)
            self.makeIntJson(eval(file.readlines()[0].replace("\n","")))
            file.close()

    def makeIntJson(self, json) :
        for key,value in json.iteritems() :
            self.json[int(key)] = value
        self.maxRunInJson = max(self.json.keys())
        self.runs = self.json.keys()

    def inJson(self) :
        run = self.source["run"]
        if self.acceptFutureRuns and run>self.maxRunInJson : return True
        if not (run in self.runs) : return False
        lumiRanges = self.json[run]
        ls = self.source["lumiSection"]
        for lumiRange in lumiRanges :
            if (ls>=lumiRange[0] and ls<=lumiRange[1]) : return True
        return False

    def update(self, ignored) :
        self.value = 1.0 if self.inJson() else None
#####################################
class PtSorted(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["Pt"])
    def update(self,_) :
        pt = self.source[self.Pt]
        self.value = all(i>j for i,j in zip(pt[:-1],pt[1:]))
#####################################
class Covariance(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(['SigmaXX','SigmaXY','SigmaYY'])
    def update(self,_) :
        self.value = np.array([[self.source[self.SigmaXX],self.source[self.SigmaXY]],
                               [self.source[self.SigmaXY],self.source[self.SigmaYY]]])
#####################################
class TriDiscriminant(wrappedChain.calculable) :
    def __init__(self, fixes = ("","") , LR = None, splitLR = 0.5, LC = None, RC = None) :
        self.fixes = fixes
        self.split = splitLR
        self.LR = LR
        self.LC = LC
        self.RC = RC

    def update(self,_) :
        R  = int( self.split < self.source[self.LR] )
        self.value = ( 1 - self.source[self.RC] ) if R else ( self.source[self.LC] - 1 )
#####################################
class KarlsruheDiscriminant(wrappedChain.calculable) :
    def __init__(self, jet, met) :
        self.stash(['M3'],jet)
        self.met = met
        self.moreName = "-8.met if met<40 else m3; %s%s; %s"%(jet+(met,))
    def update(self,_) :
        met = self.source[self.met].pt()
        self.value = -8*met if met<40 else self.source[self.M3]
