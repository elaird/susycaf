from supy import wrappedChain,utils
import configuration
import math,re,os,ROOT as r
try:
    import numpy as np
except:
    pass

#####################################
class ecalDeadTowerIsBarrel(wrappedChain.calculable) :
    etaBE = configuration.detectorSpecs()["cms"]["etaBE"]
    def update(self,ignored) : self.value = map( self.isBarrel, self.source["ecalDeadTowerTrigPrimP4"] )
    def isBarrel(self, p4) : return abs(p4.eta()) < self.etaBE
#####################################
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
class metPlusIndices(wrappedChain.calculable) :
    @property
    def name(self) :
        return "%sPlus%sIndices%s"%(self.met, self.collection[0],self.collection[1])
    def __init__(self, met, collection = None) :
        self.collection = collection
        self.fixes = collection
        self.met = met
        self.stash(["Indices","P4"])
    def update(self, ignored) :
        particleIndices = self.source[self.Indices]
        particles       = self.source[self.P4]
        self.value = self.source[self.met]
        for iParticle in particleIndices:
            self.value += particles.at(iParticle)
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
#####################################
class hcalLaserEvent2012(wrappedChain.calculable):
    def __init__(self, file="AllBadHCALLaser.txt") :
        self.file = file
        self.moreName = "pass HCAL laser 2012 event filter"
        self.setup()

    def setup(self) :
        f = file("data/%s"%self.file, 'r')
        self.laserEvents = f.read()
        f.close()

    def update(self,_):
        run = self.source["run"]
        lum = self.source["lumiSection"]
        evtNum = self.source["event"]
        if run > 203742 : self.value = 1  #HCAL calibration laser was off in RunD
        else:
            evt = "%s:%s:%s"%(run,lum,evtNum)
            laserEvent = re.search(evt, self.laserEvents)
            self.value = 0 if laserEvent else 1
#####################################
class ecalLaserCalibEvent2012(wrappedChain.calculable):
    def __init__(self) :
        self.moreName = "pass Ecal laser calibration 2012 event filter"
        self.setup()

    def setup(self) :
        f_A = file("data/ecalLaserFilter_HT_Run2012A.txt", 'r')
        f_B = file("data/ecalLaserFilter_HTMHT_Run2012B.txt", 'r')
        self.laserEvents_A = f_A.read()
        self.laserEvents_B = f_B.read()
        f_A.close()
        f_B.close()

    def update(self,_):
        run = self.source["run"]
        lum = self.source["lumiSection"]
        evtNum = self.source["event"]
        evt = "%s:%s:%s"%(run,lum,evtNum)

        if run > 196600: #last run in HTMHTParked_Run2012B
            self.value = 1
        elif run > 193621: #last run in HT_Run2012A
            laserEvent = re.search(evt, self.laserEvents_B)
            self.value = 0 if laserEvent else 1
        else:
            laserEvent = re.search(evt, self.laserEvents_A)
            self.value = 0 if laserEvent else 1
#####################################
class singleIsolatedTrack(wrappedChain.calculable):
    def __init__(self, ptMin=10., etaMax=None, dzMax = 1.0, relIso=.12, muons=None, electrons=None):
        self.ptMin = ptMin
        self.etaMax = etaMax
        self.dzMax = dzMax
        self.relIso = relIso
        self.muons = "Indices".join(muons)
        self.muonsP4 = "P4".join(muons)
        self.electrons = "Indices".join(electrons)
        self.electronsP4 = "P4".join(electrons)
        self.jetsP4 = "CorrectedP4".join(("xcak5Jet","Pat"))
        self.jets = "Indices".join(("xcak5Jet","Pat"))
        self.iEvent = 0

    def update(self,_):
        self.iEvent +=1
        #may not have p4s
        pt = self.source["PFCandsPt"]
        p4s = self.source["PFCandsP4"]
        pfCharge = self.source["PFCandsChrg"]
        pfDz = self.source["PFCandsDzPV"]
        pfTrkIso = self.source["PFCandsTrkIso"]

        muons = self.source[self.muons]
        electrons = self.source[self.electrons]
        jets = self.source[self.jets]

        mP4 = self.source[self.muonsP4]
        eP4 = self.source[self.electronsP4]
        jP4 = self.source[self.jetsP4]

        mP4List = [mP4.at(x) for x in muons]
        eP4List = [eP4.at(x) for x in electrons]
        jP4List = [jP4.at(x) for x in jets]

        self.value = []

        assert len(pt)==len(p4s)
        for ip4,p4 in enumerate(p4s):
            dr =[]
            if p4.pt() < self.ptMin: continue
            if self.etaMax:
                if p4.eta() > self.etaMax: continue
            if pfCharge[ip4] == 0: continue
            if pfDz[ip4] > self.dzMax: continue
            if pfTrkIso[ip4]/p4.pt() > self.relIso: continue
            if self.matchesIn(p4, mP4List): continue
            if self.matchesIn(p4, eP4List): continue
            self.value.append(p4)

    def matchesIn(self, p4, P4List):
        dR = 0.5
        if not P4List : return False
        matches = []
        for P4 in P4List:
            if 0.5 > r.Math.VectorUtil.DeltaR(P4,p4) :
                matches.append(P4)
        return matches
