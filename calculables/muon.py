from supy import wrappedChain,utils,calculables
import ROOT as r
##############################
class NumberOfMatches(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["IsTrackerMuon"])
        self.moreName = "hard-coded to 2"
    def isFake(self) :
        return True
    def update(self,ignored) :
        self.value = [2] * self.source[self.IsTrackerMuon].size()
##############################
class NumberOfValidPixelHits(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["IsTrackerMuon"])
        self.moreName = "hard-coded to 1"
    def isFake(self) :
        return True
    def update(self,ignored) :
        self.value = [1]*self.source[self.IsTrackerMuon].size()
##############################
class ID_TOPPAG(wrappedChain.calculable) :
    '''https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Muon_Selection'''
    def __init__(self, collection = None ) :
        self.fixes = collection
        self.stash(["IsGlobalMuon",
                    "IsTrackerMuon",
                    # pt > 20
                    # |eta|<2.1
                    "GlobalTracknormalizedChi2",
                    "GlobalTracknumberOfValidTrackerHits",
                    "GlobalTracknumberOfValidHits",
                    "GlobalTrackDxy",
                    # reliso
                    # dr mu jet
                    "Vertex",
                    "NumberOfPixelLayersWithMeasurement",
                    "NumberOfMatches"
                    ])
        self.moreName = "TOP PAG l+jets recommendation sans pt,eta,reliso,dr_(mu,jet)"

    def id(self, isGlobal, isTracker, nChi2, hitsTracker, hitsTotal, dxy, vertex, pixelLayers, nStationsMatch ) :
        return ( isGlobal and
                 isTracker and
                 nChi2 < 10 and
                 hitsTracker > 10 and
                 hitsTotal > hitsTracker and
                 dxy < 0.02 and
                 abs(vertex.z() - self.pvz) < 1 and
                 pixelLayers >= 1 and
                 nStationsMatch > 1 )
    
    def update(self,_) :
        vIndices = self.source["vertexIndices"]
        self.pvz = self.source["vertexPosition"][vIndices[0]].z() if len(vIndices) else 0
        self.value = utils.hackMap( self.id,
                                    self.source[self.IsGlobalMuon],
                                    self.source[self.IsTrackerMuon],
                                    self.source[self.GlobalTracknormalizedChi2],
                                    self.source[self.GlobalTracknumberOfValidTrackerHits],
                                    self.source[self.GlobalTracknumberOfValidHits],
                                    self.source[self.GlobalTrackDxy],
                                    self.source[self.Vertex],
                                    self.source[self.NumberOfPixelLayersWithMeasurement],
                                    self.source[self.NumberOfMatches] )
##############################
class IDtight(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["IsTrackerMuon","IDGlobalMuonPromptTight","NumberOfMatches",
                    "InnerTrackNumberOfValidHits","NumberOfValidPixelHits","GlobalTrackDxy"])
        self.moreName = "implemented by hand, CMS AN-2010/211"

    def tight(self,isTrk, idGlbTight, nStationsMatch, nTrkPxHits, nPxHits, dxy) :
        return isTrk               and \
               idGlbTight          and \
               nStationsMatch >  1 and \
               nTrkPxHits     > 10 and \
               nPxHits        >  0 and \
               abs(dxy)       <  0.2#cm

    def update(self,ignored) :
        self.value = utils.hackMap(self.tight,
                         self.source[self.IsTrackerMuon],
                         self.source[self.IDGlobalMuonPromptTight],
                         self.source[self.NumberOfMatches],
                         self.source[self.InnerTrackNumberOfValidHits],
                         self.source[self.NumberOfValidPixelHits],
                         self.source[self.GlobalTrackDxy])    
##############################
class IdPog2012Tight(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["IsGlobalMuon","IsPFMuon","GlobalTracknormalizedChi2", "GlobalTracknumberOfValidMuonHits",
                    "NumberOfMatchedStations", "InnerTrackDxy", "InnerTrackDz",
                    "NumberOfValidPixelHits", "NumberOfTrackerLayersWithMeasurement"
                    ])
        self.moreName = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon"

    def update(self,ignored) :
        self.value = []

        gm = self.source[self.IsGlobalMuon]
        pf = self.source[self.IsPFMuon]
        chi2 = self.source[self.GlobalTracknormalizedChi2]
        mHits = self.source[self.GlobalTracknumberOfValidMuonHits]
        mStations = self.source[self.NumberOfMatchedStations]
        dxy = self.source[self.InnerTrackDxy]
        dz = self.source[self.InnerTrackDz]
        pHits = self.source[self.NumberOfValidPixelHits]
        tLayers = self.source[self.NumberOfTrackerLayersWithMeasurement]

        for i in range(gm.size()) :
            reqs = [gm.at(i), pf.at(i), chi2.at(i)<10.0, mHits.at(i)>0,
                    mStations.at(i)>1, dxy.at(i)<0.2, dz.at(i)<0.5,
                    pHits.at(i)>0, tLayers.at(i)>5]
            self.value.append(all(reqs))
##############################
class CombinedRelativeIso(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["TrackIso","EcalIso","HcalIso","P4"])
        self.moreName = "(trackIso + ecalIso + hcalIso) / pt_mu"

    def combinedRelativeIso(self,isoTrk,isoEcal,isoHcal,p4) :
        return (isoTrk+isoEcal+isoHcal)/p4.pt() if p4.pt() > 0.1 else 1e10

    def update(self,ignored) :
        self.value = utils.hackMap(self.combinedRelativeIso,
                         self.source[self.TrackIso],
                         self.source[self.EcalIso],
                         self.source[self.HcalIso],
                         self.source[self.P4])
##############################
class TrackRelIso(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["TrackIso","P4"])
    def trackreliso(self,iso,p4) : return iso/(iso+p4.pt()) if iso or p4.pt() else 1
    def update(self,ignored) :
        self.value = utils.hackMap(self.trackreliso, self.source[self.TrackIso],self.source[self.P4] )
##############################
class HcalRelIso(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["HcalIso","P4"])
    def hcalreliso(self,iso,p4) : return iso/(iso+p4.pt()) if iso or p4.pt() else 1
    def update(self,ignored) :
        self.value = utils.hackMap(self.hcalreliso, self.source[self.HcalIso],self.source[self.P4] )
##############################
class IndicesOther(calculables.IndicesOther) :
    def __init__(self, collection = None) :
        super(IndicesOther, self).__init__(collection)
        self.moreName = "pass ptMin; fail id"
##############################
class IndicesNonIso(calculables.IndicesOther) :
    def __init__(self, collection = None) :
        super(IndicesNonIso, self).__init__(collection)
        self.indicesOther = "%sIndicesNonIso%s"%collection
        self.moreName = "pass ptMin & id; fail iso"
##############################
class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, isoMax = None, requireIsGlobal = True , ID = "IDtight", ISO = "CombinedRelativeIso", absEtaMax = 1000) :
        self.fixes = collection
        self.requireIsGlobal = requireIsGlobal
        self.stash(["IndicesNonIso","IndicesOther","P4","IsGlobalMuon","IsTrackerMuon","IsPFMuon"])
        self.ID = ID.join(collection)
        self.looseID = ("IdPog2012Loose").join(collection)
        self.ISO = ISO.join(collection)

        self.ptMin = ptMin
        self.absEtaMax = absEtaMax
        self.isoMax = isoMax
        self.moreName = "%s; pt>%.1f GeV; %s<%.2f"%( ID, ptMin, ISO, isoMax )

    def update(self,ignored) :
        self.value = []
        nonIso = self.source[self.IndicesNonIso]
        other  = self.source[self.IndicesOther]
        p4s    = self.source[self.P4]
        id  = self.source[self.ID]
        looseId = self.source[self.looseID]
        iso = self.source[self.ISO]
        isGlobal = self.source[self.IsGlobalMuon]
        isPF =  self.source[self.IsPFMuon]
        isTrkMu =  self.source[self.IsTrackerMuon]
        for i in range(p4s.size()) :
            p4 = p4s.at(i)
            if p4.pt() < self.ptMin : continue
            if self.requireIsGlobal and not isGlobal.at(i) : continue
            if self.absEtaMax < abs(p4.eta()) : continue
            if id[i] :
                if iso[i] < self.isoMax :
                    self.value.append(i)
                else: nonIso.append(i)
            else: other.append(i)
##############################
class IndicesAnyIso(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","IndicesNonIso"])
        self.moreName = "sorted(Indices+IndicesNonIso)"

    def update(self, ignored) :
        self.value = sorted(self.source[self.Indices]+self.source[self.IndicesNonIso])
##############################
class IndicesAnyIsoIsoOrder(wrappedChain.calculable) :
    def __init__(self, collection = None, key = None) :
        self.fixes = collection
        self.stash(["IndicesAnyIso"])
        self.key = ("%s"+key+"%s")%collection
        self.moreName = "sorted(IndicesAnyIso,key=%s)"%key

    def update(self, ignored) :
        self.value = sorted(self.source[self.IndicesAnyIso],key = self.source[self.key].__getitem__)
##############################
class LeadingPt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","P4"])

    def update(self,ignored) :
        indices = self.source[self.Indices]
        self.value = 0 if not indices else self.source[self.P4][indices[0]].pt()
class LeadingPtAny(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["IndicesAnyIso","P4"])

    def update(self,ignored) :
        indices = self.source[self.IndicesAnyIso]
        self.value = 0 if not indices else self.source[self.P4][indices[0]].pt()
class LeadingIsoAny(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, iso = None) :
        self.fixes = collection
        self.stash(["IndicesAnyIsoIsoOrder","P4"])
        self.iso = ("%s"+iso+"%s")%self.fixes
        self.ptMin = ptMin
        self.moreName = "%s of most isolated w/ pt>%.1f"%(iso,ptMin)

    def update(self,ignored) :
        p4 = self.source[self.P4]
        indices = filter(lambda i: p4[i].pt()>self.ptMin, self.source[self.IndicesAnyIsoIsoOrder])
        self.value = 10e10 if not indices else self.source[self.iso][indices[0]]
##############################
class DiMuon(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","P4"])
    def update(self,ignored) :
        indices = self.source[self.Indices]
        if len(indices)!= 2 :
            self.value = None
        else :
            p4s = self.source[self.P4]
            self.value = p4s.at(indices[0])+p4s.at(indices[1])
##############################
class DiMuonNonIso(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","P4","IndicesNonIso"])
    def update(self,ignored) :
        indices = self.source[self.Indices]
        indicesNonIso = self.source[self.IndicesNonIso]
        self.value = []
        if len(indices) + len(indicesNonIso) < 2 :
            self.value=[]
        else :
            p4s = self.source[self.P4]
            for iNonIsoMu in range(len(indicesNonIso)) :
                self.value.append(p4s.at(indices[0])+p4s.at(indicesNonIso[iNonIsoMu]))
##############################
class DiMuonOther(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["Indices","P4","IndicesOther"])
    def update(self,ignored) :
        indices = self.source[self.Indices]
        indicesOther = self.source[self.IndicesOther]
        self.value = []
        if len(indices) + len(indicesOther) < 2 :
            self.value = []
        else :
            p4s = self.source[self.P4]
            for iOtherMu in range(len(indicesOther)) :
                self.value.append(p4s.at(indices[0])+p4s.at(indicesOther[iOtherMu]))
##############################
class DiMuonMass(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["DiMuon"])
    def update(self,ignored) :
        Z = self.source[self.DiMuon]
        self.value = 0 if not Z else Z.mass()
##############################
class DiMuonNonIsoInZMass(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["DiMuonNonIso"])
    def update(self,ignored) :
        Z = self.source[self.DiMuonNonIso]
        self.value = []
        for z in Z:
            if 66.2 < z.mass() < 116.2:
                self.value.append(z.mass())
##############################
class DiMuonOtherInZMass(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["DiMuonOther"])
    def update(self,ignored) :
        Z = self.source[self.DiMuonOther]
        self.value = []
        for z in Z:
            if 66.2 < z.mass() < 116.2:
                self.value.append(z.mass())
##############################
class DiMuonPt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["DiMuon"])
    def update(self,ignored) :
        Z = self.source[self.DiMuon]
        self.value = 0 if not Z else Z.pt()
##############################
class TriggeringIndex(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(['IndicesTriggering'])
    def update(self,_) : self.value = next(iter(self.source[self.IndicesTriggering]), None)
##############################
class IndicesTriggering(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, absEtaMax = 2.1) :
        self.fixes = collection
        self.stash(["P4","IndicesAnyIso"])
        self.absEtaMax = absEtaMax
        self.ptMin = ptMin
        self.moreName = "%s in |eta|<%.1f; max pt>%.1f"%(self.IndicesAnyIso,self.absEtaMax,self.ptMin if self.ptMin else 0.)

    def update(self,ignored) :
        self.value = []
        p4 = self.source[self.P4]
        for i in self.source[self.IndicesAnyIso] :
            if p4[i].pt() < self.ptMin : break
            if abs(p4[i].eta()) < self.absEtaMax : self.value.append( i )
#####################################
class TriggeringPt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["P4","TriggeringIndex"])

    def update(self,ignored) :
        index = self.source[self.TriggeringIndex]
        self.value = 0 if index==None else self.source[self.P4][index].pt()
#####################################
class Pt(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["P4"])
    def update(self,_) :
        p4 = self.source[self.P4]
        self.value = [p4[i].pt() for i in range(len(p4))]
#####################################
class MinJetDR(wrappedChain.calculable) :
    def __init__(self, collection = None, jets = None, jetPtMin = 0) :
        self.fixes = collection
        self.ptMin = jetPtMin
        self.stash(["P4"])
        self.stash(["CorrectedP4","Indices"], jets)
        self.moreName = "min DR of mu with any jet; %s%s pt>%d"%(jets+(jetPtMin,))
        
    def update(self,_) :
        self.value = []
        mu = self.source[self.P4]
        jet = self.source[self.CorrectedP4]
        for iMu in range(len(mu)) :
            self.value.append( min([r.Math.VectorUtil.DeltaR(mu[iMu],jet[iJet]) for iJet in self.source[self.Indices] if self.ptMin < jet[iJet].pt() ]+[5]) )
##############################
class IndicesIsoLoose(wrappedChain.calculable) :
    def __init__(self,collection = None, ptMin = None, absEtaMax = None , iso = "PFIsoRel", isoMax = None) :
        self.fixes = collection
        self.stash(["IsGlobalMuon","P4"])
        self.ptMin = ptMin
        self.absEtaMax = absEtaMax
        self.iso = iso.join(collection)
        self.isoMax = isoMax
    
    def update(self,_) :
        self.value = []
        p4s = self.source[self.P4]
        isos = self.source[self.iso]
        isGlobal = self.source[self.IsGlobalMuon]
        for i in range(len(p4s)) :
            p4 = p4s[i]
            if p4.pt() < self.ptMin : break
            if ( abs(p4.eta()) < self.absEtaMax and
                 isos[i] < self.isoMax and
                 isGlobal[i]
                 ) : self.value.append(i)
##############################
class IdPog2012Loose(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["IsGlobalMuon","IsTrackerMuon","IsPFMuon",
                    ])
        self.moreName = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon"

    def update(self,ignored) :
        self.value = []

        gm = self.source[self.IsGlobalMuon]
        pf = self.source[self.IsPFMuon]
        tr = self.source[self.IsTrackerMuon]

        for i in range(gm.size()) :
            if pf.at(i):
                if gm.at(i) or tr.at(i) :
                    self.value.append(i)
            else:
                break
##############################
