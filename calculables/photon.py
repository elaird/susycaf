from supy import wrappedChain,calculables,utils
import configuration
import collections
##############################
class IndicesOther(calculables.IndicesOther) :
    def __init__(self,collection = None) :
        super(IndicesOther, self).__init__(collection)
        self.moreName = "pass ptMin; fail id/iso"
##############################
class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, flagName = None):
        self.fixes = collection
        self.stash(["IndicesOther", "P4"])
        self.ptMin = ptMin
        self.flagName = flagName
        self.moreName = "pT>=%.1f GeV; %s"% (ptMin, flagName if flagName else "")
    def update(self,ignored) :
        p4s = self.source[self.P4]
        ids = self.source[self.flagName] if self.flagName else p4s.size()*[1]
        self.value = []
        other = self.source[self.IndicesOther]
        for i in range(p4s.size()):
            if p4s.at(i).pt() < self.ptMin: continue
            elif ids[i] : self.value.append(i)
            else: other.append(i)
##############################
class leadingPt(wrappedChain.calculable) :
    @property
    def name(self) : return "%sLeadingPt%s"%self.photons

    def __init__(self, collection = ("photon","Pat")) :
        self.photons = collection
        
    def update(self,ignored) :
        indices = self.source["%sIndices%s"%self.photons]
        self.value = self.source["%sP4%s"%self.photons].at(indices[0]).pt() if len(indices) else None
####################################
class photonID(wrappedChain.calculable) :
    @property
    def name(self) : return self.idName
    
    # following https://twiki.cern.ch/twiki/bin/viewauth/CMS/PhotonID
    def __init__(self, collection = None, level = None) :
        self.cs = collection
        self.idName = "%sID%s%s" % (self.cs[0],level,self.cs[1])
        self.p4Name = "%sP4%s" % self.cs

        for var in ["EcalRecHitEtConeDR04", "HcalTowSumEtConeDR04",
                    "HadronicOverEm", "TrkSumPtHollowConeDR04", "SigmaIetaIeta","HasPixelSeed"] :
            setattr(self,var, ("%s"+var+"%s")%self.cs)

        jei = {}; jeiLower = {}; tbhi = {}; tbhiLower = {}; hcti = {}; hctiLower = {};
        hoe = {}; shhBarrel = {}; shhEndcap = {}; ptVar = {}; etaBE = {}; moreName = {}
        for l in ["EmFromTwiki","LooseFromTwiki","TightFromTwiki",
                  "AnalysisNote_10_268","EGM_10_006_Loose","EGM_10_006_Tight","TrkIsoSideBand","TrkIsoRelaxed","IsoSideBand","IsoRelaxed","NoIsoReq"] :
            jei [l]      = (4.2,  0.0060); jeiLower[l]  = None
            tbhi[l]      = (2.2,  0.0025); tbhiLower[l] = None
            hcti[l]      = None          ; hctiLower[l] = None
            hoe [l]      = (0.05, 0.0000)
            shhBarrel[l] = None
            shhEndcap[l] = None
            ptVar[l]     = "pt"
            etaBE[l]     = configuration.detectorSpecs()["cms"]["etaBE"]
            moreName[l]  = "PhotonID twiki, 2010-10-14, %s"%("is"+l.replace("FromTwiki",""))

        hcti     ["LooseFromTwiki"] = (3.5, 0.001)

        hcti     ["TightFromTwiki"] = (2.0, 0.001)
        shhBarrel["TightFromTwiki"] = (0.013, 0.0)
        shhEndcap["TightFromTwiki"] = (0.030, 0.0)

        jei      ["AnalysisNote_10_268"] = (4.2, 0.003)
        tbhi     ["AnalysisNote_10_268"] = (2.2, 0.001)
        hcti     ["AnalysisNote_10_268"] = (2.0, 0.001)
        ptVar    ["AnalysisNote_10_268"] = "Et"
        moreName ["AnalysisNote_10_268"] = "from CMS AN 10-268"

        jei      ["EGM_10_006_Tight"] = (2.4,  0.0)
        tbhi     ["EGM_10_006_Tight"] = (1.0,  0.0)
        hcti     ["EGM_10_006_Tight"] = (0.9,  0.0)
        hoe      ["EGM_10_006_Tight"] = (0.03, 0.0)
        shhBarrel["EGM_10_006_Tight"] = (0.01, 0.0)
        shhEndcap["EGM_10_006_Tight"] = (0.028,0.0)
        moreName ["EGM_10_006_Tight"] = "EGM 10-006 tight"

        jei      ["EGM_10_006_Loose"] = (4.2,  0.0)
        tbhi     ["EGM_10_006_Loose"] = (2.2,  0.0)
        hcti     ["EGM_10_006_Loose"] = (2.0,  0.0)
        hoe      ["EGM_10_006_Loose"] = (0.05, 0.0)
        shhBarrel["EGM_10_006_Loose"] = (0.01, 0.0)
        shhEndcap["EGM_10_006_Loose"] = (0.03, 0.0)
        moreName ["EGM_10_006_Loose"] = "EGM 10-006 loose"

        hcti     ["TrkIsoRelaxed"] = (10.0, 0.001)
        moreName ["TrkIsoRelaxed"] = "relaxed trkIso"

        hcti     ["TrkIsoSideBand"] = (10.0, 0.001)
        hctiLower["TrkIsoSideBand"] = ( 2.0, 0.001)
        moreName ["TrkIsoSideBand"] = "side-band of trkIso"

        jei      ["NoIsoReq"] = (100.0, 0.0)
        tbhi     ["NoIsoReq"] = (100.0, 0.0)
        hcti     ["NoIsoReq"] = (100.0, 0.0)
        moreName ["NoIsoReq"] = "relaxed trkIso [ ,100]; hcalIso[ ,100]; ecalIso[ ,100]"

        jei      ["IsoRelaxed"] = (8.2,  0.0060)
        tbhi     ["IsoRelaxed"] = (6.2,  0.0025)
        hcti     ["IsoRelaxed"] = (10.0, 0.001)
        moreName ["IsoRelaxed"] = "relaxed trkIso [ ,10]; hcalIso[ ,6]; ecalIso[ ,8]"

        jei      ["IsoSideBand"] = (8.2,  0.0060)
        jeiLower ["IsoSideBand"] = jei ["TightFromTwiki"]
        tbhi     ["IsoSideBand"] = (6.2,  0.0025)
        tbhiLower["IsoSideBand"] = tbhi["TightFromTwiki"]
        hcti     ["IsoSideBand"] = (10.0, 0.001)
        hctiLower["IsoSideBand"] = hcti["TightFromTwiki"]
        moreName ["IsoSideBand"] = "side-band of trkIso [2,10]; hcalIso[2,6]; ecalIso[4,8]"

        for item in ["jei","jeiLower",
                     "tbhi","tbhiLower",
                     "hcti","hctiLower",
                     "hoe","shhBarrel","shhEndcap","ptVar","etaBE","moreName"] :
            setattr(self,item,eval(item)[level])

    def update(self,ignored) :
        self.value = utils.hackMap(self.passId, 
                         self.source[self.p4Name],
                         self.source[self.EcalRecHitEtConeDR04],
                         self.source[self.HcalTowSumEtConeDR04],
                         self.source[self.HadronicOverEm],
                         self.source[self.TrkSumPtHollowConeDR04],
                         self.source[self.SigmaIetaIeta],
                         self.source[self.HasPixelSeed],
                         )

    def passId(self, p4, jei, tbhi, hoe, hcti, shh, hasPixelSeed) :
        def pass1(item, value, alsoLower = False) :
            member = getattr(self,item)
            if member!=None and value > (member[0] + pt*member[1]) :
                return False
            if alsoLower :
                memberLower = getattr(self,item+"Lower")
                if memberLower!=None and value < (memberLower[0] + pt*memberLower[1]) : return False
            return True
        
        pt = getattr(p4,self.ptVar)()

        if not pass1("jei",  eval("jei"),  alsoLower = True)  : return False
        if not pass1("tbhi", eval("tbhi"), alsoLower = True)  : return False
        if not pass1("hcti", eval("hcti"), alsoLower = True)  : return False
        if not pass1("hoe",  eval("hoe"),  alsoLower = False) : return False

        shhVar = self.shhBarrel if abs(p4.eta())<self.etaBE else self.shhEndcap
        if shhVar!=None and shh  > (shhVar[0] + pt*shhVar[1]) : return False
        
        if hasPixelSeed : return False
        return True
####################################    
class HcalTowSumEtConeDR04(wrappedChain.calculable) :
    @property
    def name(self) : return "%sHcalTowSumEtConeDR04%s"%self.collection
    def __init__(self, collection = None) :
        self.collection = collection
        self.var1 = "%sHcalDepth1TowSumEtConeDR04%s"%collection
        self.var2 = "%sHcalDepth2TowSumEtConeDR04%s"%collection
        self.moreName = "depth 1 + depth 2"
    def update(self, ignored) :
        size = len(self.source[self.var1])
        self.value = [self.source[self.var1].at(i)+self.source[self.var2].at(i) for i in range(size)]
####################################    
class SeedTime(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["P4"])
    def isFake(self) : return True
    def update(self, ignored) :
        self.value = [-100.0]*len(self.source[self.P4])
####################################    
class WrappedSeedTime(wrappedChain.calculable) :
    def __init__(self, collection) :
        self.fixes = collection
        self.stash(["P4","SeedTime"])
    def update(self, ignored) :
        nPhot  = len(self.source[self.P4])
        nSeeds = len(self.source[self.SeedTime])
        if nPhot!=nSeeds :
            self.value = [-100.0]*len(self.source[self.P4])
        else :
            self.value = self.source[self.SeedTime]
####################################
class CombinedIsoDR03RhoCorrected(wrappedChain.calculable) :
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["EcalRecHitEtConeDR03", "HcalTowSumEtConeDR03", "TrkSumPtHollowConeDR03"])

    def update(self, _) :
        self.value = []
        e = self.source[self.EcalRecHitEtConeDR03]
        h = self.source[self.HcalTowSumEtConeDR03]
        t = self.source[self.TrkSumPtHollowConeDR03]
        rho = self.source["rho"]
        for i in range(e.size()) :
            self.value.append(e.at(i) + h.at(i) + t.at(i) - rho*(0.093+0.0281))
####################################
class IDRA3(wrappedChain.calculable) :
    #https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=188427
    def __init__(self, collection = None) :
        self.fixes = collection
        self.stash(["CombinedIsoDR03RhoCorrected", "HadronicOverEm", "SigmaIetaIeta", "HasPixelSeed", "R9"])

    def update(self, _) :
        self.value = []
        hoe = self.source[self.HadronicOverEm]
        shh = self.source[self.SigmaIetaIeta]
        pix = self.source[self.HasPixelSeed]
        r9  = self.source[self.R9]
        iso = self.source[self.CombinedIsoDR03RhoCorrected]

        for i in range(hoe.size()) :
            self.value.append(all([not pix.at(i), hoe.at(i)<0.05, shh.at(i)<0.011, r9.at(i)<1.0, iso[i]<6.0]))
####################################
class photonIDEmFromTwiki(photonID) :
    def __init__(self, collection = None) :
        super(photonIDEmFromTwiki,self).__init__(collection,"EmFromTwiki")
####################################
class photonIDLooseFromTwiki(photonID) :
    def __init__(self, collection = None) :
        super(photonIDLooseFromTwiki,self).__init__(collection,"LooseFromTwiki")
####################################
class photonIDTightFromTwiki(photonID) :
    def __init__(self, collection = None) :
        super(photonIDTightFromTwiki,self).__init__(collection,"TightFromTwiki")
####################################
class photonIDTrkIsoSideBand(photonID) :
    def __init__(self, collection = None) :
        super(photonIDTrkIsoSideBand,self).__init__(collection,"TrkIsoSideBand")
####################################
class photonIDNoIsoReq(photonID) :
    def __init__(self, collection = None) :
        super(photonIDNoIsoReq,self).__init__(collection,"NoIsoReq")
####################################
class photonIDIsoRelaxed(photonID) :
    def __init__(self, collection = None) :
        super(photonIDIsoRelaxed,self).__init__(collection,"IsoRelaxed")
####################################
class photonIDIsoSideBand(photonID) :
    def __init__(self, collection = None) :
        super(photonIDIsoSideBand,self).__init__(collection,"IsoSideBand")
####################################
class photonIDTrkIsoRelaxed(photonID) :
    def __init__(self, collection = None) :
        super(photonIDTrkIsoRelaxed,self).__init__(collection,"TrkIsoRelaxed")
####################################
class photonIDAnalysisNote_10_268(photonID) :
    def __init__(self, collection = None) :
        super(photonIDAnalysisNote_10_268,self).__init__(collection,"AnalysisNote_10_268")
####################################
class photonIDEGM_10_006_Tight(photonID) :
    def __init__(self, collection = None) :
        super(photonIDEGM_10_006_Tight,self).__init__(collection,"EGM_10_006_Tight")
####################################
class photonIDEGM_10_006_Loose(photonID) :
    def __init__(self, collection = None) :
        super(photonIDEGM_10_006_Loose,self).__init__(collection,"EGM_10_006_Loose")
####################################
class photonWeight(wrappedChain.calculable) :
    def __init__(self, var = None, weightSet = "November 2011") :
        self.var = var

        self.weight = collections.defaultdict(int)
        if weightSet=="June 2011" :
            self.weight.update({0:1.00002, 1:0.246676, 2:0.784571, 3:1.34209, 4:1.49489, 5:1.45355, 6:1.10074, 7:0.90572, 8:0.595763, 9:0.514093, 10:0.469262, 11:0.451208, 12:0.59306, 13:0.459517, 14:0.467984, 15:0, 16:0, 17:0, 18:0, 19:0})#determined on 2011-06-14 requiring a photon

        if weightSet=="November 2011" :
            self.weight.update({0:1.00004, 1:0.0872381, 2:0.441893, 3:0.764065, 4:1.08425, 5:1.2493, 6:1.34958, 7:1.44752, 8:1.47106, 9:1.60048, 10:1.75302, 11:1.79496, 12:2.28168, 13:2.33737, 14:2.2478, 15:3.22463, 16:4.26286, 17:4.71759, 18:7.22803, 19:24.0934})#determined on 2011-11-29 requiring a photon

        if weightSet=="ZM" :
            lst = [ 0.0, 0.737935, 6.20516, 4.81946, 7.17781, 5.23649, 5.87606, 6.50951, 4.63189, 4.57838, 3.53186, 2.93797, 2.11116, 1.64937, 1.40657, 0.936336, 0.748872, 0.591696, 0.407115, 0.317994, 0.243221, 0.187518, 0.115325, 0.151001, 0.0843434, 0.0516356, 0.049791, 0.0272625, 0.0189111, 0.00867613, 0.00481305, 0, 0.0184357, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0, 0, 0, 0, 0 ]
            for i,w in enumerate(lst) :
                self.weight[i] = w

    def update(self, ignored) :
        self.value = self.weight[len(self.source[self.var])]
####################################
class photonWeightChoppedToOne(wrappedChain.calculable) :
    def __init__(self, var = None) :
        self.weight = {0:0.999993, 1:7.7039, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0, 18:0, 19:0}#determined on 2011-05-05
        self.var = var
    def update(self, ignored) :
        self.value = self.weight[len(self.source[self.var])]
####################################
