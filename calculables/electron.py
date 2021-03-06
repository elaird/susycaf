from supy import wrappedChain,utils
import configuration
from muon import IndicesOther,IndicesNonIso,IndicesAnyIso,IndicesAnyIsoIsoOrder,LeadingPt
##############################
barrelEtaMax = configuration.detectorSpecs()["cms"]["barrelEtaMax"]
endcapEtaMin = configuration.detectorSpecs()["cms"]["endcapEtaMin"]
##############################
class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, simpleEleID = None, useCombinedIso = True, flag2012 = "") :
        self.fixes = collection
        self.stash(["IndicesOther","IndicesNonIso","P4"])
        self.ptMin = ptMin
        self.flag2012 = flag2012
        if self.flag2012 :
            self.eID = ("%sId"+self.flag2012+"%s")% self.fixes
            self.moreName = "pt>%.1f; (2012 %s)" % (self.ptMin, self.flag2012)
        else :
            isotype = "c" if useCombinedIso else "rel"
            self.eID = ("%sID"+simpleEleID+"%s")% self.fixes
            self.eIso = ("%s"+isotype+"Iso"+simpleEleID+"%s") % self.fixes
            self.moreName = "pt>%.1f; simple%s; %sIso" % (ptMin, simpleEleID, isotype)

    def update(self,ignored) :
        self.value = []
        other = self.source[self.IndicesOther]
        nonIso = self.source[self.IndicesNonIso]
        p4s   = self.source[self.P4]
        eID   = self.source[self.eID]
        eIso  = self.source[self.eIso] if not self.flag2012 else [True]*p4s.size()
        for i in range(p4s.size()) :
            if p4s.at(i).pt() < self.ptMin : break
            elif eID[i]:
                if eIso[i]: self.value.append(i)
                else : nonIso.append(i)
            else: other.append(i)
##############################
class ID(wrappedChain.calculable) :
    def id(self, hoe, dphi, deta, sieie, missingHits, dist, DcotTh, scEta) :
        absEta = abs(scEta)
        notConv = missingHits <= self.missingHits and (self.dist==None or dist < self.dist) and (self.DcotTh==None or DcotTh < self.DcotTh)
        notConv|= not self.rejectConversions
        return notConv and hoe < self.Bhoe and dphi < self.Bdphi and deta < self.Bdeta and sieie < self.Bsieie if absEta < barrelEtaMax else \
               notConv and hoe < self.Ehoe and dphi < self.Edphi and deta < self.Edeta and sieie < self.Esieie if absEta > endcapEtaMin else \
               False

    def update(self,ignored) :
        self.value = utils.hackMap(self.id, 
                         self.source[self.HcalOverEcal],
                         self.source[self.DeltaPhiSuperClusterTrackAtVtx],
                         self.source[self.DeltaEtaSuperClusterTrackAtVtx],
                         self.source[self.SigmaIetaIeta],
                         self.source[self.ConversionMissingHits],
                         self.source[self.ConversionDist],
                         self.source[self.ConversionDCot],
                         self.source[self.ESuperClusterEta])

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID
    def __init__(self, collection, eff, rejectConversions = True) :
        self.fixes = collection
        self.stash(["HcalOverEcal","DeltaPhiSuperClusterTrackAtVtx","DeltaEtaSuperClusterTrackAtVtx","SigmaIetaIeta",
                    "ConversionDist", "ConversionDCot", "ConversionMissingHits","ESuperClusterEta"])
        i = ["95","90","85","80","70","60"].index(eff)
        self.rejectConversions = rejectConversions
        
        self.Bsieie = [0.01,    0.01,    0.01,    0.01,    0.01,    0.01  ][i]
        self.Bdphi  = [0.80,    0.80,    0.06,    0.06,    0.03,    0.025 ][i]
        self.Bdeta  = [0.007,   0.007,   0.006,   0.004,   0.004,   0.004 ][i]
        self.Bhoe   = [0.15,    0.12,    0.04,    0.04,    0.025,   0.025 ][i]
        self.Esieie = [0.03,    0.03,    0.03,    0.03,    0.03,    0.03  ][i]
        self.Edphi  = [0.7,     0.7,     0.04,    0.03,    0.02,    0.02  ][i]
        self.Edeta = 10.0 #[0.01,   0.009,  0.007,  0.007,  0.005,  0.005 ][i] # note on twiki not to apply these
        self.Ehoe   = [0.07,    0.05,    0.025,   0.025,   0.025,   0.025 ][i]
        #conversion rejection
        self.missingHits = [1,  1,       1,       0,       0,       0     ][i]
        self.dist   = [None,    0.02,    0.02,    0.02,    0.02,    0.02  ][i]
        self.DcotTh = [None,    0.02,    0.02,    0.02,    0.02,    0.02  ][i]
        
class ID95NoConversionRejection(ID) :
    def __init__(self,collection) : super(ID95NoConversionRejection,self).__init__(collection,"95", rejectConversions = False)
class ID95(ID) :
    def __init__(self,collection) : super(ID95,self).__init__(collection,"95")
class ID90(ID) :
    def __init__(self,collection) : super(ID90,self).__init__(collection,"90")
class ID85(ID) :
    def __init__(self,collection) : super(ID85,self).__init__(collection,"85")
class ID80(ID) :
    def __init__(self,collection) : super(ID80,self).__init__(collection,"80")
class ID70(ID) :
    def __init__(self,collection) : super(ID70,self).__init__(collection,"70")
class ID60(ID) :
    def __init__(self,collection) : super(ID60,self).__init__(collection,"60")
##############################
class Iso(wrappedChain.calculable) :
    def __init__(self, collection, eff, combined) :
        self.fixes = collection
        self.stash(["IsoCombined","P4","TrackIsoRel","EcalIsoRel","HcalIsoRel","ESuperClusterEta"])
        self.combined = combined
        self.isoName = collection[0]+("c" if combined else "rel")+"Iso"+eff+collection[1]
        i = ["95","90","85","80","70","60"].index(eff)
        self.Bc    = [0.15,   0.10,   0.09,   0.07,   0.04,   0.03 ][i]
        self.Btrk  = [0.15,   0.12,   0.09,   0.09,   0.05,   0.04 ][i]
        self.Becal = [2.00,   0.09,   0.08,   0.07,   0.06,   0.04 ][i]
        self.Bhcal = [0.12,   0.10,   0.10,   0.10,   0.03,   0.03 ][i]
        self.Ec    = [0.1,    0.07,   0.06,   0.06,   0.03,   0.02 ][i]
        self.Etrk  = [0.08,   0.05,   0.05,   0.04,   0.025,  0.025][i]
        self.Eecal = [0.06,   0.06,   0.05,   0.05,   0.025,  0.02 ][i]
        self.Ehcal = [0.05,   0.03,   0.025,  0.025,  0.02,   0.02 ][i]
        
    def cIso(self, iso, scEta) :
        absEta = abs(scEta)
        if absEta < barrelEtaMax: return iso < self.Bc
        if absEta > endcapEtaMin: return iso < self.Ec
        return None

    def relIso(self, trk, ecal, hcal, scEta) :
        absEta = abs(scEta)
        if absEta < barrelEtaMax : return trk < self.Btrk and ecal < self.Becal and hcal < self.Bhcal
        if absEta > endcapEtaMin : return trk < self.Etrk and ecal < self.Eecal and hcal < self.Ehcal
        return None

    def update(self,ignored) :
        self.value = utils.hackMap(self.cIso,
                         self.source[self.IsoCombined],
                         self.source[self.ESuperClusterEta]) if self.combined else \
                     utils.hackMap(self.relIso,
                         self.source[self.TrackIsoRel],
                         self.source[self.EcalIsoRel],
                         self.source[self.HcalIsoRel],
                         self.source[self.ESuperClusterEta])
class cIso95NoConversionRejection(Iso) :
    def __init__(self,collection) : super(cIso95NoConversionRejection,self).__init__(collection,"95", True)
class cIso95(Iso) :
    def __init__(self,collection) : super(cIso95,self).__init__(collection,"95", True)
class cIso90(Iso) :
    def __init__(self,collection) : super(cIso90,self).__init__(collection,"90", True)
class cIso85(Iso) :
    def __init__(self,collection) : super(cIso85,self).__init__(collection,"85", True)
class cIso80(Iso) :
    def __init__(self,collection) : super(cIso80,self).__init__(collection,"80", True)
class cIso70(Iso) :
    def __init__(self,collection) : super(cIso70,self).__init__(collection,"70", True)
class cIso60(Iso) :
    def __init__(self,collection) : super(cIso60,self).__init__(collection,"60", True)
class relIso95(Iso) :
    def __init__(self,collection) : super(relIso95,self).__init__(collection,"95", False)
class relIso90(Iso) :
    def __init__(self,collection) : super(relIso90,self).__init__(collection,"90", False)
class relIso85(Iso) :
    def __init__(self,collection) : super(relIso85,self).__init__(collection,"85", False)
class relIso80(Iso) :
    def __init__(self,collection) : super(relIso80,self).__init__(collection,"80", False)
class relIso70(Iso) :
    def __init__(self,collection) : super(relIso70,self).__init__(collection,"70", False)
class relIso60(Iso) :
    def __init__(self,collection) : super(relIso60,self).__init__(collection,"60", False)
##############################
class IsoCombined(wrappedChain.calculable) :
    def __init__(self,collection = None) :
        self.fixes = collection
        self.stash(["Dr03TkSumPt","Dr03EcalRecHitSumEt","Dr03HcalTowerSumEt","P4","ESuperClusterEta"])

    def update(self,ignored) :
        self.value = utils.hackMap(self.combinedIso,
                                   self.source[self.Dr03TkSumPt],
                                   self.source[self.Dr03EcalRecHitSumEt],
                                   self.source[self.Dr03HcalTowerSumEt],
                                   self.source[self.P4],
                                   self.source[self.ESuperClusterEta])

    def combinedIso(self,trk,ecal,hcal,p4,scEta) :
        absEta = abs(scEta)
        return (trk + max(0.,ecal-1) + hcal) / p4.pt() if absEta < barrelEtaMax else \
               (trk + ecal + hcal) / p4.pt() if absEta > endcapEtaMin else \
               None
##############################
class IsoCombinedAdjusted(wrappedChain.calculable) :
    '''IsoCombined, scaled only for endcap by constant M, such that barrelCIso == M * endcapCIso.

    Default values of {barrelCIso:1,endcapCIso:1} provide
    IsoCombinedAdjusted values identical to those of IsoCombined.
    '''

    def __init__(self,collection = None, barrelCIso = 1, endcapCIso = 1) :
        self.fixes = collection
        self.stash(['IsoCombined','ESuperClusterEta'])
        self.endcapScaling = barrelCIso / endcapCIso

    def update(self,_) :
        isos = self.source[self.IsoCombined]
        scetas = self.source[self.ESuperClusterEta]
        self.value = [ iso if abs(scetas.at(i)) < barrelEtaMax else
                       iso * self.endcapScaling if iso!=None else None
                       for i,iso in enumerate(isos)]
##############################
class IsoRel(wrappedChain.calculable) :
    def __init__(self, collection, isoSource) :
        self.fixes = collection
        self.stash(["P4"])
        self.isoSource = ("%s"+isoSource+"%s") % collection
    def relIso(self, iso, p4) : return iso/p4.pt()
    def update(self,ignored) : self.value = utils.hackMap(self.relIso, self.source[self.isoSource], self.source[self.P4])
    
class TrackIsoRel(IsoRel) :
    def __init__(self, collection = None) : super(TrackIsoRel,self).__init__(collection, "Dr03TkSumPt")
class EcalIsoRel(IsoRel) :
    def __init__(self, collection = None) : super(EcalIsoRel,self).__init__(collection, "Dr03EcalRecHitSumEt")
class HcalIsoRel(IsoRel) :
    def __init__(self, collection = None) : super(HcalIsoRel,self).__init__(collection, "Dr03HcalTowerSumEt")
##############################
class ConversionMissingHits(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["HcalOverEcal"])
        self.moreName = "hard-coded to 0"
    def isFake(self) :
        return True
    def update(self, ignored) :
        self.value = [0]*len(self.source[self.HcalOverEcal])
##############################
class Indices_TopPAG(wrappedChain.calculable) :
    @property
    def name(self) : return "Indices".join(self.fixes)

    def __init__(self, collection = None, ptMin = None, absEtaMax = None, id = "ID70") :
        self.fixes = collection
        self.ptMin = ptMin
        self.absEtaMax = absEtaMax
        self.id = id.join(collection)
        self.stash(["P4","GsfTrackDxyBS","GsfTrackVertexz",])
        
    def update(self,_) :
        p4 = self.source[self.P4]
        vz = self.source[self.GsfTrackVertexz]
        dxyBS = self.source[self.GsfTrackDxyBS]

        vIndices = self.source["vertexIndices"]
        pvz = self.source["vertexPosition"][vIndices[0]].z() if len(vIndices) else 0
        
        self.value = [i for i,id in enumerate(self.source[self.id])
                      if ( id and
                           self.ptMin < p4.at(i).pt() and
                           abs(p4.at(i).eta()) < self.absEtaMax and
                           abs(vz.at(i) - pvz) < 1.0 and
                           abs(dxyBS.at(i)) < 0.02 )]
##############################
class IndicesIsoLoose(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, absEtaMax = None, iso = "PFIsoRel", isoMax = None) :
        self.fixes = collection
        self.stash(["P4"])
        self.ptMin = ptMin
        self.absEtaMax = absEtaMax
        self.iso = iso.join(collection)
        self.isoMax = isoMax

    def update(self,_) :
        self.value = []
        p4s = self.source[self.P4]
        isos = self.source[self.iso]
        for i in range(len(p4s)) :
            p4 = p4s.at(i)
            if p4.pt() < self.ptMin : break
            if ( isos[i] != None and
                 abs(p4.eta()) < self.absEtaMax and
                 isos[i] < self.isoMax
                 ) : self.value.append(i)                 
##############################
