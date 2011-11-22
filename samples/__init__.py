for module in [
    "mc",
    "jetMet",
    "ht",
    "muon",
    "electron",
    "photon",
    "signalSkim",
    "wPol",
    "doubleMu",
    ] : exec("from __%s__ import %s"%s)
