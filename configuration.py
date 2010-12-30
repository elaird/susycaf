def calculablesFiles() :
    return [("Gen","gen"), ("Jet","jet"), ("Muon","muon"), ("Electron","electron"), ("Photon","photon"), ("Other","other"), ("XClean","xclean"), ("Compatibility","compat")]

def samplesFiles() :
    return ["MC", "JetMET", "Muon", "Photon", "SignalSkim"]

def stepsFiles() :
    return ["Other", "Jet", "Trigger", "Photon", "Print", "Gen", "Xclean", "Displayer"]

def stepsToDisableForData() :
    return ["genMotherHistogrammer", "photonPurityPlots", "photonEfficiencyPlots"]

def histogramsToDisableForData() :
    return ["genpthat"]

def stepsToDisableForMc() :
    return ["hltFilter", "hltFilterList", "lowestUnPrescaledTrigger", "lowestUnPrescaledTriggerHistogrammer", "hbheNoiseFilter", "bxFilter", "physicsDeclared", "techBitFilter"]

def histogramsToDisableForMc() :
    return []

def sourceFiles() :
    return ["pragmas.h","helpers.C"]

def batchScripts(hostName) :
    d = {"lx05.hep.ph.ic.ac.uk":("icSub.sh","icJob.sh"),
         "lx06.hep.ph.ic.ac.uk":("icSub.sh","icJob.sh"),
         }
    assert hostName in d,"hostname %s not recognized"%hostName
    return d[hostName]
