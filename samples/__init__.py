for module in [
    "mc",
    "mcOld",
    "jetmet",
    "ht",
    "muon",
    "electron",
    "photon",
    "signalSkim",
    "wpol",
    "mumu",
    "top",
    "ewk",
    ] : exec("from __%s__ import %s"%(module,module))
