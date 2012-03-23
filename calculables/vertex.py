import ROOT as r
from supy import wrappedChain,utils,calculables
#####################################
class ID(wrappedChain.calculable) :
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2010Recipes#Good_Vertex_selection
    def __init__(self, minNdof = 5.0, maxAbsZ = 24.0, maxD0 = 2.0 ) :
        self.fixes = ("vertex","")
        for item in ["minNdof","maxAbsZ","maxD0"] : setattr(self,item,eval(item))
        self.moreName = "!fake; nd>=%.1f; |z|<=%.1f cm; d0<=%.1f cm" % (minNdof,maxAbsZ,maxD0)

    def id(self, isFake, ndof, pos) :
        return (not isFake) and \
               ndof >= self.minNdof and \
               abs(pos.Z()) <= self.maxAbsZ and \
               abs(pos.Rho()) <= self.maxD0

    def update(self,ignored) :
        self.value = utils.hackMap(self.id, self.source["vertexIsFake"],self.source["vertexNdof"],self.source["vertexPosition"])
#####################################
class Indices(wrappedChain.calculable) :
    def __init__(self, sumPtMin = None) :
        self.fixes = ("vertex","")
        self.sumPtMin = sumPtMin
        self.moreName = ""
        if self.sumPtMin!=None :
            self.moreName += "sumPt >=%.1f"%sumPtMin
        self.moreName += "; pass ID"
        
    def update(self,ignored) :
        sumPt = self.source["vertexSumPt"]
        sumPt2 = self.source["vertexSumPt2"]
        id = self.source["vertexID"]
        self.value = []
        other = self.source["vertexIndicesOther"]
        for i in range(len(id)) :
            if self.sumPtMin!=None and sumPt.at(i) < self.sumPtMin : continue
            elif id[i] : self.value.append(i)
            else : other.append(i)
        self.value.sort( key = sumPt2.__getitem__, reverse = True )
#####################################
class IndicesOther(calculables.IndicesOther) :
    def __init__(self) :
        self.fixes = ("vertex","")
        super(IndicesOther, self).__init__(self.fixes)
        self.moreName = "pass sumPtMin; fail ID"
#####################################
class SumPt(wrappedChain.calculable) :
    def __init__(self) :
        self.fixes = ("vertex","")
        self.sumPts = r.std.vector('double')()
        for i in range(100) :
            self.sumPts.push_back(-100.0)
            
    def update(self, ignored) :
        self.value = self.sumPts
#####################################
class SumP3(wrappedChain.calculable) :
    def __init__(self) :
        self.fixes = ("vertex","")
        self.sumP3s = r.std.vector(type(utils.LorentzV()))()

    def update(self, ignored) :
        self.value = self.sumP3s
#####################################
class nVertex(wrappedChain.calculable) :
    def update(self,ignored) : self.value = len(self.source["vertexIndices"])
#####################################
class vertex0Ntracks(wrappedChain.calculable) :
    def update(self,_) :
        ntrk = self.source["vertexNtrks"]
        indices = self.source["vertexIndices"]
        self.value = ntrk[indices[0]] if len(indices) else 0
#####################################
class vertexDzSeparation(wrappedChain.calculable) :
    def update(self,_) :
        p = self.source["vertexPosition"]
        posZ = [p[i].z() for i in self.source["vertexIndices"]]
        self.value = min([29] + [abs(posZ[0]-z) for z in posZ[1:]])
#####################################
class vertexTrackedDz(wrappedChain.calculable) :
    def update(self,_) :
        p = self.source["vertexPosition"]
        ntrk = self.source["vertexNtrks"]
        indices = self.source["vertexIndices"]
        posZ0 = p[indices[0]].z()
        self.value = sum([abs(posZ0-z) for z in posZ[1:]]) / sum(ntrk[i] for i in indices[1:] )
#####################################
class vertexTrackPurity(wrappedChain.calculable) :
    def update(self,_) :
        p = self.source["vertexPosition"]
        ntrk = self.source["vertexNtrks"]
        indices = self.source["vertexIndices"]
        posZ0 = p[indices[0]].z()
        self.value = float( ntrk[indices[0]] ) / sum( ntrk[i] for i in indices if abs(posZ0-p[i].z()) < 5 )
