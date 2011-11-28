from supy import analysisStep
#####################################
class runLsEvent(analysisStep) :
    def __init__(self, fileName) :
        self.moreName = "run:ls:event in %s"%fileName
        file = open(fileName)
        self.tuples = [ eval("(%s,%s,%s)"%tuple(line.replace("\n","").split(":"))) for line in file]
        file.close()

    def select (self,eventVars) :
        return (eventVars["run"], eventVars["lumiSection"], eventVars["event"]) in self.tuples
#####################################
class run(analysisStep) :
    def __init__(self,runs,acceptRatherThanReject) :
        self.runs = runs
        self.accept = acceptRatherThanReject

        self.moreName = "run%s in list %s" % ( ("" if self.accept else " not"),str(runs) )
        
    def select (self,eventVars) :
        return not ((eventVars["run"] in self.runs) ^ self.accept)
#####################################
class hbheNoise(analysisStep) :
    def select (self,eventVars) :
        return eventVars["hbheNoiseFilterResult"]
#####################################
class monster(analysisStep) :
    def __init__(self,maxNumTracks=10,minGoodTrackFraction=0.25) :
        self.maxNumTracks=maxNumTracks
        self.minGoodTrackFraction=minGoodTrackFraction

        self.moreName = "<=%d tracks or >%.2f good fraction" % (maxNumTracks, minGoodTrackFraction)

    def select (self,eventVars) :
        nTracks = sum(map(eventVars.__getitem__, ["tracksNEtaLT0p9AllTracks",
                                                  "tracksNEta0p9to1p5AllTracks",
                                                  "tracksNEtaGT1p5AllTracks"]))
        nGoodTracks = sum(map(eventVars.__getitem__, ["tracksNEtaLT0p9HighPurityTracks",
                                                      "tracksNEta0p9to1p5HighPurityTracks",
                                                      "tracksNEtaGT1p5HighPurityTracks"]))
        return (nTracks <= self.maxNumTracks or nGoodTracks > self.minGoodTrackFraction*nTracks)
#####################################
class DeltaRGreaterSelector(analysisStep) :

    def __init__(self, jets = None, particles = None, minDeltaR = None, particleIndex = None):
        self.particleIndex = particleIndex
        self.minDeltaR = minDeltaR
        self.particleIndicesName = "%sIndices%s"%particles

        self.moreName = "%s%s; DR(%s%s[i[%d]], jet) > %.2f"%(jets[0], jets[1], particles[0], particles[1], particleIndex, minDeltaR)
        self.minDeltaRVar = "%s%sMinDeltaRToJet%s%s"%(particles[0], particles[1], jets[0], jets[1])
        
    def select (self,eventVars) :
        indices = eventVars[self.particleIndicesName]
        if len(indices) <= self.particleIndex : return False
        index = indices[self.particleIndex]
        return eventVars[self.minDeltaRVar][index]>self.minDeltaR
#####################################
