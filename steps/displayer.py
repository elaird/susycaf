import ROOT as r
import configuration
import supy

try:
    import pdgLookup
    pdgLookupExists = True
except ImportError:
    pdgLookupExists = False


def jetsRaw(jets):
    return (jets[0].replace("xc",""), jets[1])

class displayer(supy.steps.displayer):
    def __init__(self,
                 # objects
                 jets=None,
                 jetIndices=[],
                 met=None, muons=None, electrons=None,
                 photons=None, taus=None, recHits=None, recHitPtThreshold=-100,
                 jetsOtherAlgo=None, metOtherAlgo=None, recHitsOtherAlgo=None,
                 # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
                 bThresholds={"JetProbabilityBJetTags":         {"L":0.275,
                                                                 "M":0.545,
                                                                 "T":0.790},
                              "CombinedSecondaryVertexBJetTags":{"L":0.244,
                                                                 "M":0.679,
                                                                 "T":0.898},
                              },
                 genMet="genmetP4True",
                 genParticleIndices=[],
                 triggersToPrint=[],
                 flagsToPrint=[#"logErrorTooManyClusters",
                               #"logErrorTooManySeeds",
                               #"beamHaloCSCLooseHaloId",
                               #"beamHaloCSCTightHaloId",
                               #"beamHaloEcalLooseHaloId",
                               #"beamHaloEcalTightHaloId",
                               #"beamHaloGlobalLooseHaloId",
                               #"beamHaloGlobalTightHaloId",
                               #"beamHaloHcalLooseHaloId",
                               #"beamHaloHcalTightHaloId"
                               ],

                 # options
                 scale=200.0,
                 prettyMode=False, tipToTail=False,
                 printExtraText=True, printGen=False,
                 rename={},

                 # RA1
                 ra1Mode=True, ra1CutBits=True, j2Factor=None,
                 mhtOverMetName="", showAlphaTMet=True,
                 deltaPhiStarExtraName="", deltaPhiStarCut=None, deltaPhiStarDR=None,
                 ):

        self.moreName = "(see below)"

        for item in ["jets", "jetIndices",
                     "met", "muons", "electrons", "photons", "taus",
                     "recHits", "recHitPtThreshold", "jetsOtherAlgo",
                     "metOtherAlgo", "recHitsOtherAlgo", "bThresholds",
                     "genMet", "genParticleIndices",
                     "triggersToPrint", "flagsToPrint",

                     "scale",
                     "prettyMode", "tipToTail",
                     "printExtraText", "printGen",

                     "ra1Mode", "ra1CutBits", "j2Factor",
                     "mhtOverMetName", "showAlphaTMet",
                     "deltaPhiStarExtraName", "deltaPhiStarCut", "deltaPhiStarDR",
                     ]:
            setattr(self, item, eval(item))

        if len(self.flagsToPrint) > 3:
            print "WARNING: More than three flags specified in the displayer."+\
                  "  The list will run off the page."

        self.etaBE = configuration.detectorSpecs()["cms"]["etaBE"]
        self.subdetectors = {None: []}
        self.recHitCollections = {None: []}
        for rh in [self.recHits, self.recHitsOtherAlgo]:
            if not rh:
                continue
            self.recHitCollections[rh] = configuration.detectorSpecs()["cms"]["%sRecHitCollections" % rh]
            self.subdetectors[rh] = configuration.detectorSpecs()["cms"]["%sSubdetectors" % rh]

        if self.jets:
            self.genJets = ("gen%sGenJets" % jetsRaw(self.jets)[0][:3], "")
            self.deltaHtName = "%sDeltaPseudoJetEt%s" % self.jets
            self.deltaPhiStarName = "".join([self.jets[0], "DeltaPhiStar",
                                             self.deltaPhiStarExtraName,
                                             self.jets[1]])
            radius = 0.7 if "ak7Jet" in self.jets[0] else 0.5
            self.jetIndices = [(self.genJets, "", r.kBlack, 2, radius),
                               (self.jets, "%sIndices%s" % self.jets, r.kBlue, 1, radius),
                               (self.jets, "%sIndicesOther%s" % self.jets, r.kOrange+3, 1, radius),
                               (self.jets, "%sIndicesIgnored%s" % self.jets, r.kCyan, 1, radius),
                               ] + self.jetIndices
            #jpt = (self.jets[0].replace("xc","")+"JPT", self.jets[1])
            #pf =  (self.jets[0].replace("xc","")+"PF", self.jets[1])
            #self.jetIndices.append((jpt, "%sIndices%s" % jpt, 896, 1))
            #self.jetIndices.append((pf, "%sIndices%s" % pf, 38, 1))
        else:
            self.genJets = None

        self.prettyReName = {
            "genak5GenJetsP4": "AK5 Gen Jets",
            "xcak5JetIndicesPat": "AK5 Calo Jets",
            "xcak5JetIndicesOtherPat": "AK5 Calo Jets (fail eta or ID)",
            "xcak5JetIndicesIgnoredPat": "AK5 Calo Jets (fail pT or xc)",
            "xcak5JetIndicesBtagged2Pat": "AK5 Calo Jets (b-tagged)",

            "ak5JetPFIndicesPat": "AK5 PF Jets",
            "ak5JetPFIndicesOtherPat": "AK5 PF Jets (fail ID or eta)",
            "ak5JetPFIndicesIgnoredPat": "AK5 PF Jets (fail pT or xc)",
            "ak5JetPFIndicesBtagged2Pat": "AK5 PF Jets (b-tagged)",

            "MHT (xcak5JetPat)": "MHT",
            "MHT (xcak5JetPFPat)": "MHT",
            "metP4AK5TypeII": "MET (Calo Type II)",
            "metP4PF": "MET (PF)",
            "metP4TypeIPF": "MET (PF Type I)",
            "muonP4Pat": "muons",
            "electronP4Pat": "electrons",
            "photonP4Pat": "photons",
            }
        self.prettyReName.update(rename)

        self.titleSizeFactor = 1.0
        
        self.legendDict = {}
        self.legendList = []

        self.ellipse = r.TEllipse()
        self.ellipse.SetFillStyle(0)

        self.deadBox = r.TBox()
        self.deadBox.SetFillColor(r.kMagenta)
        self.deadBox.SetLineColor(r.kMagenta)

        self.coldBox = r.TBox()
        self.coldBox.SetFillColor(r.kOrange+7)
        self.coldBox.SetLineColor(r.kOrange+7)

        self.hcalBox = r.TBox()
        self.hcalBox.SetFillColor(r.kGreen)
        self.hcalBox.SetLineColor(r.kGreen)
        
        self.line = r.TLine()
        self.arrow = r.TArrow()
        self.text = r.TText()
        self.latex = r.TLatex()

        self.alphaFuncs=[
            self.makeAlphaTFunc(0.55,r.kBlack),
            self.makeAlphaTFunc(0.50,r.kOrange+3),
            self.makeAlphaTFunc(0.45,r.kOrange+7)
            ]

    def prepareText(self, params, coords) :
        self.text.SetTextSize(params["size"])
        self.text.SetTextFont(params["font"])
        self.text.SetTextColor(params["color"])
        self.textSlope = params["slope"]

        self.textX = coords["x"]
        self.textY = coords["y"]
        self.textCounter = 0

    def printText(self, message) :
        self.text.DrawText(self.textX, self.textY - self.textCounter * self.textSlope, message)
        self.textCounter += 1

    def printEvent(self, eventVars, params, coords) :
        self.prepareText(params, coords)
        for message in ["Run %6d / Ls %6d / Event %10d"%(eventVars["run"], eventVars["lumiSection"], eventVars["event"])
                        #"PtHat(GeV) %#5.1f"%eventVars["genpthat"] if not eventVars["isRealData"] else "",
                        ] :
            if message : self.printText(message)
        for item in self.triggersToPrint :
            self.printText("%s"%(item if eventVars["triggered"][item] else ""))
        
    def printVertices(self, eventVars, params, coords, nMax):
        if "vertexIndices" not in eventVars:
            return
        self.prepareText(params, coords)
        self.printText("Vertices")
        self.printText("ID   Z(cm)%s"%(" sumPt(GeV)" if not self.prettyMode else ""))
        self.printText("----------%s"%("-----------" if not self.prettyMode else ""))

        nVertices = eventVars["vertexNdof"].size()
        for i in range(nVertices) :
            if nMax<=i :
                self.printText("[%d more not listed]"%(nVertices-nMax))
                break
            
            out = "%2s  %6.2f"%("G " if i in eventVars["vertexIndices"] else "  ", eventVars["vertexPosition"].at(i).z())
            if not self.prettyMode : out += " %5.0f"%eventVars["vertexSumPt"].at(i)
            self.printText(out)

    def printPhotons(self, eventVars, params, coords, photons, nMax) :
        self.prepareText(params, coords)
        p4Vector = eventVars["%sP4%s"        %photons]
        loose    = eventVars["%sIDLooseFromTwiki%s"%photons]
        tight    = eventVars["%sIDTightFromTwiki%s"%photons]
            
        self.printText(self.renamedDesc(photons[0]+photons[1]))
        self.printText("ID   pT  eta  phi")
        self.printText("-----------------")

        nPhotons = p4Vector.size()
        for iPhoton in range(nPhotons) :
            if nMax<=iPhoton :
                self.printText("[%d more not listed]"%(nPhotons-nMax))
                break
            photon=p4Vector[iPhoton]

            outString = "%1s%1s"% ("L" if loose[iPhoton] else " ", "T" if tight[iPhoton] else " ")
            outString+="%5.0f %4.1f %4.1f"%(photon.pt(), photon.eta(), photon.phi())
            self.printText(outString)

    def printElectrons(self, eventVars, params, coords, electrons, nMax) :
        self.prepareText(params, coords)
        p4Vector = eventVars["%sP4%s"        %electrons]
        cIso = eventVars["%sIsoCombined%s"%electrons]
        ninetyFive = eventVars["%sID95%s"%electrons]
     
        self.printText(self.renamedDesc(electrons[0]+electrons[1]))
        self.printText("ID   pT  eta  phi")#  cIso")
        self.printText("-----------------")#------")

        nElectrons = p4Vector.size()
        for iElectron in range(nElectrons) :
            if nMax<=iElectron :
                self.printText("[%d more not listed]"%(nElectrons-nMax))
                break
            electron=p4Vector[iElectron]

            outString = "%2s"%("95" if ninetyFive[iElectron] else "  ")
            outString+="%5.0f %4.1f %4.1f"%(electron.pt(), electron.eta(), electron.phi())
            #outString+=" %5.2f"%cIso[iElectron] if cIso[iElectron]!=None else " %5s"%"-"
            self.printText(outString)

    def printMuons(self, eventVars, params, coords, muons, nMax) :
        self.prepareText(params, coords)
        p4Vector = eventVars["%sP4%s"     %muons]
        tight    = eventVars["%sIdPog2012Tight%s"%muons]
        iso      = eventVars["%sPfIsolationR04DeltaBCorrected%s"%muons]
        gl       = eventVars["%sIsGlobalMuon%s"%muons]
        pf       = eventVars["%sIsPFMuon%s"%muons]
        
        self.printText(self.renamedDesc(muons[0]+muons[1]))
        self.printText(" ID  pT  eta  phi pfIso")
        self.printText("-----------------------")

        nMuons = p4Vector.size()
        for iMuon in range(nMuons) :
            if nMax<=iMuon :
                self.printText("[%d more not listed]"%(nMuons-nMax))
                break
            muon=p4Vector[iMuon]

            outString = "%s%s%s"%("G" if gl[iMuon] else " ", "P" if pf[iMuon] else " ", "T" if tight[iMuon] else " ")
            outString+= "%4.0f %4.1f %4.1f"%(muon.pt(), muon.eta(), muon.phi())
            outString+= " %5.2f"%(iso[iMuon]) if iso[iMuon]<100.0 else ">100".rjust(6)

            self.printText(outString)

    def printRecHits(self, eventVars, params, coords, recHits, nMax) :
        self.prepareText(params, coords)

        self.printText(self.renamedDesc("%s %sRecHits"%("severe" if recHits=="Calo" else "cleaned", recHits)))
        self.printText("  det  pT  eta  phi%s"%(" sl" if recHits=="Calo" else ""))
        self.printText("-------------------%s"%("---" if recHits=="Calo" else ""))
        self.printText("SumPt%4.0f"%eventVars["%sRecHitSumPt"%recHits])
        p4 = eventVars["%sRecHitSumP4"%recHits]
        self.printText("SumP4%4.0f %s %4.1f"%(p4.pt(),
                                              "%4.1f"%p4.eta() if abs(p4.eta())<10.0 else ">10.",
                                              p4.phi()))

        hits = []
        for detector in self.subdetectors[recHits] :
            for collectionName in self.recHitCollections[recHits] :
                p4Var = "rechit%s%s%s%s"%(collectionName, recHits, "P4",            detector)
                slVar = "rechit%s%s%s%s"%(collectionName, recHits, "SeverityLevel", detector)
                for iHit in range(len(eventVars[p4Var])) :
                    hit = eventVars[p4Var].at(iHit)
                    l = [hit.pt(), hit.eta(), hit.phi(), detector]
                    if recHits=="Calo" : l.append(eventVars[slVar].at(iHit))
                    hits.append( tuple(l) )
        
        for iHit,hit in enumerate(reversed(sorted(hits))) :
            if nMax<=iHit :
                self.printText("[%d more not listed]"%(len(hits)-nMax))
                break
            outString = "%5s"%hit[3]
            outString+="%4.0f %4.1f %4.1f"%(hit[0], hit[1], hit[2])
            if recHits=="Calo" : outString +=" %2d"%hit[4]
            self.printText(outString)
        
    def printJets(self, eventVars, params, coords, jets, nMax) :
        def bCategory(key = "", value = None) :
            if key not in self.bThresholds :
                return " "
            dct = self.bThresholds[key]
            for cat in ["T","M","L"] :
                if value>dct[cat] : return cat
            return " "

        self.prepareText(params, coords)
        jets2 = jetsRaw(jets)
        isPf = "PF" in jets[0]
        
        p4Vector         = eventVars['%sCorrectedP4%s'%jets]
        corrVector       = eventVars['%sCorrFactor%s'%jets2]
        csv              = eventVars['%sCombinedSecondaryVertexBJetTags%s'%jets2]
        jp               = eventVars['%sJetProbabilityBJetTags%s'%jets2]
        bIndices         = eventVars['%sIndicesBtagged2%s'%jets] if '%sIndicesBtagged2%s'%jets in eventVars else []
        
        if not isPf :
            jetEmfVector  = eventVars['%sEmEnergyFraction%s'%jets2]
            jetFHpdVector = eventVars['%sJetIDFHPD%s'       %jets2]
            jetFRbxVector = eventVars['%sJetIDFRBX%s'       %jets2]
            jetN90Vector  = eventVars['%sJetIDN90Hits%s'    %jets2]
            
            loose = eventVars["%sJetIDloose%s"%jets2]
            tight = eventVars["%sJetIDtight%s"%jets2]
            
        else :
            chf = eventVars["%sFchargedHad%s"%jets2]
            nhf = eventVars["%sFneutralHad%s"%jets2]

            cef = eventVars["%sFchargedEm%s"%jets2]
            nef = eventVars["%sFneutralEm%s"%jets2]

            cm  = eventVars["%sNcharged%s"%jets2]
            nm  = eventVars["%sNneutral%s"%jets2]
            
            loose = eventVars["%sPFJetIDloose%s"%jets2]
            tight = eventVars["%sPFJetIDtight%s"%jets2]
            
        self.printText(self.renamedDesc(jets[0]+jets[1]))
        left = "ID    pT  eta  phi"
        center = "   EMF  fHPD  fRBX N90" if not isPf else "   CHF  NHF  CEF  NEF CM"
        #right = " corr   CSV"
        right = "   CSV     JP"
        self.printText(left+center+right)
        self.printText("-"*(len(left)+len(center)+len(right)))

        nJets = p4Vector.size()
        for iJet in range(nJets) :
            if nMax<=iJet :
                self.printText("[%d more not listed]"%(nJets-nMax))
                break
            jet=p4Vector[iJet]

            outString = "%1s%1s%1s"% ("L" if loose[iJet] else " ", "T" if tight[iJet] else " ", "b" if iJet in bIndices else " ")
            outString+="%5.0f %4.1f %4.1f"%(jet.pt(), jet.eta(), jet.phi())

            if not isPf :
                outString+=" %5.2f %5.2f %5.2f %3d"%(jetEmfVector.at(iJet), jetFHpdVector.at(iJet), jetFRbxVector.at(iJet), jetN90Vector.at(iJet))
            else :
                outString+=" %5.3f %4.2f %4.2f %4.2f%3d"%(chf.at(iJet), nhf.at(iJet), cef.at(iJet), nef.at(iJet), cm.at(iJet))

            csvValue = csv.at(iJet)
            #outString += " %4.2f %7.3f"%(corrVector.at(iJet), csv.at(iJet))

            if csvValue==-10. :
                outString += " %4.0f"%csvValue
            else :
                outString += " %4.2f"%csvValue
            outString += bCategory("CombinedSecondaryVertexBJetTags", csvValue)

            jpValue = jp.at(iJet)
            outString += " %5.2f"%jpValue
            outString += bCategory("JetProbabilityBJetTags", jpValue)
            self.printText(outString)

    def printGenJets(self, eventVars, params, coords, nMax) :
        p4Key = "%sP4%s" % self.genJets
        if p4Key not in eventVars:
            return
        p4Vector = eventVars[p4Key]
            
        self.prepareText(params, coords)
        self.printText(self.renamedDesc(p4Key))
        self.printText("   pT  eta  phi")
        self.printText("---------------")

        nJets = p4Vector.size()
        for iJet in range(nJets):
            if nMax <= iJet:
                self.printText("[%d more not listed]"%(nJets-nMax))
                break
            jet = p4Vector.at(iJet)
            self.printText("%5.0f %4.1f %4.1f"%(jet.pt(), jet.eta(), jet.phi()))

    def printGenParticles(self, eventVars, params, coords, nMax) :
        self.prepareText(params, coords)
        p4s    = eventVars["genP4"]
        status = eventVars["genStatus"]
        ids    = eventVars["genPdgId"]
        
        self.printText("Status 1 Gen Particles")
        self.printText("  name  pdgId   pT  eta  phi")
        self.printText("----------------------------")

        particles = reversed(sorted([(i, p4s[i]) for i in range(p4s.size())], key = lambda x:x[1].pt()))
        nStatus1 = sum([status[i]==1 for i in range(status.size())])
        iPrint = 0
        for iParticle,p4 in particles :
            if status.at(iParticle)!=1 : continue
            if nMax<=iPrint :
                self.printText("[%d more not listed]"%(nStatus1-nMax))
                break
            pdgId = ids.at(iParticle)
            self.printText("%6s %6d%5.0f %4.1f %4.1f"%(pdgLookup.pdgid_to_name(pdgId) if pdgLookupExists else "", pdgId, p4.pt(), p4.eta(), p4.phi()))
            iPrint += 1
        return
    
    def printKinematicVariables(self, eventVars, params, coords, jets, jets2) :
        self.prepareText(params, coords)
        
        def go(j) :
            dps = eventVars[self.deltaPhiStarName]
            l = [eventVars["%sHtBin%s"%j],
                 eventVars["%s%s%s"  %(j[0], "SumEt",        j[1])],
                 eventVars["%s%s%s"  %(j[0], "SumP4",        j[1])].pt() if eventVars["%s%s%s"%(j[0], "SumP4",  j[1])] else 0,
                 eventVars["%s%s%s"  %(j[0], "AlphaTEt",     j[1])],
                 dps[0][0] if dps else -1.0,
                 ]
            for i in range(len(l)) :
                if l[i]==None : l[i] = -1.0
            self.printText("%14s %4.0f %4.0f %4.0f %6.3f %5.2f"%tuple([self.renamedDesc(j[0]+j[1])]+l))

        self.printText("jet collection  bin   HT  MHT alphaT Dphi*")
        self.printText("------------------------------------------")
        
        go(jets)
        if jets2!=None :
            go(jets2)
        
    def passBit(self, var) :
        return " p" if var else " f"

    def printCutBits(self, eventVars, params, coords, jets, jets2, met, met2) :
        self.prepareText(params, coords)

        def go(j, m, i) :
            J2 = None if len(eventVars["%sIndices%s"%j])<2 else eventVars['%sCorrectedP4%s'%j].at(eventVars["%sIndices%s"%j][1]).pt()
            HT = eventVars["%sSumEt%s"%j]
            aT = eventVars["%sAlphaTEt%s"%j]
            MM = eventVars[self.mhtOverMetName]
            dedr = eventVars["%sDeadEcalDR%s%s"%(j[0], self.deltaPhiStarExtraName, j[1])]
            DE = (not dedr) or dedr[0]>self.deltaPhiStarDR

            htBin = None
            if eventVars["%sHtBin%s"%j] : htBin = eventVars["%sHtBin%s"%j]
            elif eventVars["%sFixedHtBin%s"%j] : htBin = eventVars["%sFixedHtBin%s"%j]

            j2Bit = J2!=None and htBin!=None and J2 > self.j2Factor*htBin
            htBit = HT!=None and htBin!=None and HT > htBin
            atBit = aT!=None and aT > 0.550
            deBit = DE!=None and DE
            mmBit = MM!=None and MM < 1.250
            
            self.printText("%14s  %s %s %s %s %s  %s"%(self.renamedDesc(j[0]+j[1]),
                                                       self.passBit(j2Bit),
                                                       self.passBit(htBit),
                                                       self.passBit(atBit),
                                                       self.passBit(deBit),
                                                       self.passBit(mmBit),
                                                       "candidate (%g)"%htBin if all([j2Bit, htBit, atBit, deBit, mmBit]) else "",
                                                       )
                           )
            if self.prettyMode and all and not i :
                self.text.SetTextSize(1.5*params["size"])
                self.text.SetTextFont(params["font"])
                self.text.SetTextColor(r.kBlue)
                #self.text.DrawText(0.1, 0.1, "passes final selection")
                self.text.SetTextSize(params["size"])
                self.text.SetTextFont(params["font"])
                self.text.SetTextColor(params["color"])

        self.printText("jet collection  J2 HT aT DE MM")
        self.printText("------------------------------")

        go(jets, met, 0)
        if jets2!=None and met2!=None :
            go(jets2, met2, 1)

    def printFlags(self, eventVars, params, coords, flags) :
        self.prepareText(params, coords)
        for f in flags:
            if eventVars[f] : self.printText(f)
        
    def drawSkeleton(self, coords, color) :
        r.gPad.AbsCoordinates(False)
        
        self.ellipse.SetLineColor(color)
        self.ellipse.SetLineWidth(1)
        self.ellipse.SetLineStyle(1)
        self.ellipse.DrawEllipse(coords["x0"], coords["y0"], coords["radius"], coords["radius"], 0.0, 360.0, 0.0, "")

        self.line.SetLineColor(color)
        self.line.DrawLine(coords["x0"]-coords["radius"], coords["y0"]                 , coords["x0"]+coords["radius"], coords["y0"]                 )
        self.line.DrawLine(coords["x0"]                 , coords["y0"]-coords["radius"], coords["x0"]                 , coords["y0"]+coords["radius"])

    def drawScale(self, color, size, scale, point) :
        self.latex.SetTextSize(size)
        self.latex.SetTextColor(color)
        self.latex.DrawLatex(point["x"], point["y"],"radius = "+str(scale)+" GeV p_{T}")

    def drawP4(self,
               coords=None,
               p4=None,
               lineColor=None,
               lineWidth=1,
               lineStyle=1,
               arrowSize=1.0,
               p4Initial=None):

        c = coords
        x0 = c["x0"]+p4Initial.px()*c["radius"]/c["scale"] if p4Initial else c["x0"]
        y0 = c["y0"]+p4Initial.py()*c["radius"]/c["scale"] if p4Initial else c["y0"]
        x1 = x0+p4.px()*c["radius"]/c["scale"]
        y1 = y0+p4.py()*c["radius"]/c["scale"]

        #self.line.SetLineColor(lineColor)
        #self.line.SetLineWidth(lineWidth)
        #self.line.DrawLine(x0,y0,x1,y1)

        self.arrow.SetLineColor(lineColor)
        self.arrow.SetLineWidth(lineWidth)
        self.arrow.SetLineStyle(lineStyle)
        self.arrow.SetArrowSize(arrowSize)
        self.arrow.SetFillColor(lineColor)
        self.arrow.DrawArrow(x0,y0,x1,y1)
        
    def drawCircle(self, p4=None, lineColor=None, lineWidth=1, lineStyle=1, circleRadius=None):
        self.ellipse.SetLineColor(lineColor)
        self.ellipse.SetLineWidth(lineWidth)
        self.ellipse.SetLineStyle(lineStyle)
        self.ellipse.DrawEllipse(p4.eta(), p4.phi(), circleRadius, circleRadius, 0.0, 360.0, 0.0, "")

    def renamedDesc(self, desc):
        if not self.prettyMode:
            return desc
        desc = desc.replace("IndicesStatus3", "Status3")
        return desc if desc not in self.prettyReName else self.prettyReName[desc]

    def legendFunc(self, lineColor=None, lineStyle=1, name="", desc=""):
        if name not in self.legendDict:
            self.legendDict[name] = True
            self.legendList.append((lineColor, lineStyle, self.renamedDesc(desc), "l"))

    def drawGenParticles(self, eventVars=None, indices="",
                         coords=None, lineColor=None,
                         lineWidth=1, lineStyle=1,
                         arrowSize=-1.0, circleRadius=None):

        self.legendFunc(lineColor=lineColor,
                        lineStyle=lineStyle,
                        name=indices,
                        desc=indices)

        for iParticle in eventVars[indices]:
            particle = eventVars["genP4"].at(iParticle)
            if circleRadius is None:
                self.drawP4(coords=coords,
                            p4=particle,
                            lineColor=lineColor,
                            lineWidth=lineWidth,
                            arrowSize=arrowSize)
            else :
                self.drawCircle(p4=particle,
                                lineColor=lineColor,
                                lineWidth=lineWidth,
                                circleRadius=circleRadius)


    def drawJets(self, eventVars=None, jets=None, indices="",
                 coords=None, lineColor=None, lineWidth=1, lineStyle=1,
                 arrowSize=-1.0, circleRadius=None):

        p4Key = '%s%sP4%s' % (jets[0], "" if "gen" in jets[0] else 'Corrected', jets[1])
        if p4Key not in eventVars:
            return

        self.legendFunc(lineColor=lineColor,
                        lineStyle=lineStyle,
                        name=indices,
                        desc=indices if indices else p4Key)
        p4s = eventVars[p4Key]
        if self.tipToTail:
            phiOrder = eventVars["%sIndicesPhi%s" % jets]
            partials = eventVars["%sPartialSumP4%s" % jets]
            mean = supy.utils.partialSumP4Centroid(partials)
            for i in range(len(phiOrder)) :
                self.drawP4(coords=coords,
                            p4=p4s.at(phiOrder[i]),
                            lineColor=lineColor,
                            lineWidth=lineWidth,
                            lineStyle=lineStyle,
                            arrowSize=0.3*arrowSize,
                            p4Initial=partials[i]-mean)
            return

        for iJet in eventVars[indices] if indices else range(len(p4s)):
            p4 = p4s.at(iJet)
            if circleRadius:
                self.drawCircle(p4=p4,
                                lineColor=lineColor,
                                lineStyle=lineStyle,
                                circleRadius=circleRadius)
            else:
                self.drawP4(coords=coords,
                            p4=p4,
                            lineColor=lineColor,
                            lineWidth=lineWidth,
                            lineStyle=lineStyle,
                            arrowSize=arrowSize)


    def drawMht(self, eventVars, coords, color, lineWidth, arrowSize) :
        self.legendFunc(lineColor=color, name="%smht%s"%self.jets, desc="MHT (%s%s)"%self.jets)

        sump4 = eventVars["%sSumP4%s"%self.jets]
        if self.tipToTail :
            phiOrder = eventVars["%sIndicesPhi%s"%self.jets]
            partials = eventVars["%sPartialSumP4%s"%self.jets]
            mean = eventVars["%sPartialSumP4Centroid%s"%self.jets]
            if sump4 : self.drawP4(coords, -sump4,color,lineWidth,
                                   arrowSize=arrowSize, p4Initial=partials[-1]-mean)
            return
        if sump4 : self.drawP4(coords, -sump4, color, lineWidth, arrowSize=arrowSize)
            
    def drawHt(self, eventVars, coords, color, lineWidth, arrowSize) :
        self.legendFunc(lineColor=color,
                        name="%sht%s" % self.jets,
                        desc="H_{T} (%s%s)" % self.jets)

        ht = eventVars["%sSumPt%s"%self.jets]
            
        y = coords["y0"]-coords["radius"]-0.05
        l = ht*coords["radius"]/coords["scale"]
        self.line.SetLineColor(color)
        self.line.DrawLine(coords["x0"]-l/2.0, y, coords["x0"]+l/2.0, y)
        
    def drawNJetDeltaHt(self, eventVars, coords, color, lineWidth, arrowSize) :
        self.legendFunc(lineColor=color,
                        name="%sdeltaHt%s" % self.jets,
                        desc="#DeltaH_{T} (%s%s)" % self.jets)

        y = coords["y0"]-coords["radius"]-0.03
        l = eventVars[self.deltaHtName]*coords["radius"]/coords["scale"]
        self.line.SetLineColor(color)
        self.line.DrawLine(coords["x0"]-l/2.0, y, coords["x0"]+l/2.0, y)

    def drawGenMet(self, eventVars, coords, color, lineWidth, arrowSize) :
        if self.genMet==None : return
        self.legendFunc(lineColor=color,
                        name="genMet",
                        desc="GEN MET (%s)" % self.genMet)
        self.drawP4(coords, eventVars[self.genMet], color, lineWidth, arrowSize=arrowSize)
            
    def drawCleanedRecHits(self, eventVars, coords, color, lineWidth, arrowSize) :
        self.legendFunc(lineColor=color,
                        name="cleanedRecHits%s" % self.recHits,
                        desc="cleaned RecHits (%s)" % self.recHits)

        for detector in self.subdetectors :
            for collectionName in self.recHitCollections :
                varName = "rechit%s%sP4%s"%(collectionName, self.recHits, detector)
                for iHit in range(len(eventVars[varName])) :
                    hit = eventVars[varName].at(iHit)
                    if hit.pt()<self.recHitPtThreshold : continue
                    self.drawP4(coords, hit, color, lineWidth, arrowSize=arrowSize)

    def drawCleanedRecHitSumP4(self, eventVars, coords, color, lineWidth, arrowSize) :
        self.legendFunc(lineColor=color,
                        name="%sRecHitSumP4" % self.recHits,
                        desc="severe %sRechits SumP4" % self.recHits)
        sump4 = eventVars["%sRecHitSumP4"%self.recHits]
        if sump4 : self.drawP4(coords, sump4, color, lineWidth, arrowSize=arrowSize)
            
    def makeAlphaTFunc(self,alphaTValue,color) :
        alphaTFunc=r.TF1("#alpha_{T} = %#4.2g"%alphaTValue,
                         "1.0-2.0*("+str(alphaTValue)+")*sqrt(1.0-x*x)",
                         0.0,1.0)
        alphaTFunc.SetLineColor(color)
        alphaTFunc.SetLineWidth(1)
        alphaTFunc.SetNpx(300)
        return alphaTFunc

    def etaPhiPad(self, eventVars, corners):
        pad = r.TPad("etaPhiPad", "etaPhiPad",
                     corners["x1"], corners["y1"],
                     corners["x2"], corners["y2"])
        pad.cd()
        pad.SetTickx()
        pad.SetTicky()

        etaPhiPlot = r.TH2D("etaPhi", ";#eta;#phi;",
                            1, -r.TMath.Pi(), r.TMath.Pi(),
                            1, -r.TMath.Pi(), r.TMath.Pi())
        etaPhiPlot.SetStats(False)
        etaPhiPlot.Draw()

        self.line.SetLineColor(r.kBlack)
        self.line.DrawLine(-self.etaBE, etaPhiPlot.GetYaxis().GetXmin(),
                           -self.etaBE, etaPhiPlot.GetYaxis().GetXmax())
        self.line.DrawLine( self.etaBE, etaPhiPlot.GetYaxis().GetXmin(),
                            self.etaBE, etaPhiPlot.GetYaxis().GetXmax())
        return pad, etaPhiPlot

    def ra1EtaPhi(self, eventVars, etaPhiPad):
        def drawEcalBox(fourVector, nBadXtals, maxStatus) :
            value = (0.087/2) * nBadXtals / 25
            args = (fourVector.eta()-value, fourVector.phi()-value, fourVector.eta()+value, fourVector.phi()+value)
            if maxStatus==14 :
                self.deadBox.DrawBox(*args)
            else :
                self.coldBox.DrawBox(*args)
                
        def drawHcalBox(fourVector) :
            value = 0.087/2
            args = (fourVector.eta()-value, fourVector.phi()-value, fourVector.eta()+value, fourVector.phi()+value)
            self.hcalBox.DrawBox(*args)

        etaPhiPad.cd()

        #draw dead ECAL regions
        nRegions = eventVars["ecalDeadTowerTrigPrimP4"].size()
        for iRegion in range(nRegions):
            drawEcalBox(fourVector=eventVars["ecalDeadTowerTrigPrimP4"].at(iRegion),
                        nBadXtals=eventVars["ecalDeadTowerNBadXtals"].at(iRegion),
                        maxStatus=eventVars["ecalDeadTowerMaxStatus"].at(iRegion),
                        )

        #draw masked HCAL regions
        nBadHcalChannels = eventVars["hcalDeadChannelP4"].size()
        for iChannel in range(nBadHcalChannels):
            drawHcalBox(fourVector=eventVars["hcalDeadChannelP4"].at(iChannel))

        legend1 = r.TLegend(0.02, 0.9, 0.72, 1.0)
        legend1.SetFillStyle(0)
        legend1.SetBorderSize(0)
        legend1.AddEntry(self.deadBox,"dead ECAL cells","f")
        legend1.AddEntry(self.coldBox,"dead ECAL cells w/TP link","f")
        legend1.AddEntry(self.hcalBox,"masked HCAL cells","f")
        legend1.Draw()

        # FIXME: make this a calculable
        suspiciousJetIndices = []
        for dPhiStar, iJet in eventVars[self.deltaPhiStarName]:
            if dPhiStar < self.deltaPhiStarCut:
                suspiciousJetIndices.append(iJet)

        suspiciousJetColor = r.kRed
        suspiciousJetStyle = 2
        suspiciousJetLegendEntry = False
        p4s = eventVars["%sCorrectedP4%s" % self.jets]
        for iJet in suspiciousJetIndices:
            self.drawCircle(p4=p4s.at(iJet),
                            lineColor=suspiciousJetColor,
                            circleRadius=self.deltaPhiStarDR,
                            lineStyle=suspiciousJetStyle)
            suspiciousJetLegendEntry = True

        legend2 = r.TLegend(0.48, 0.933, 0.98, 1.0)
        legend2.SetFillStyle(0)
        legend2.SetBorderSize(0)
        self.ellipse.SetLineColor(suspiciousJetColor)
        self.ellipse.SetLineStyle(suspiciousJetStyle)
        if suspiciousJetLegendEntry:
            legend2.AddEntry(self.ellipse,
                             "jet with min. #Delta#phi* < %3.1f" % self.deltaPhiStarCut,
                             "l")
        legend2.Draw()

        return [legend1, legend2]


    def drawAlphaPlot(self, eventVars, color, corners):
        pad = r.TPad("alphaTpad", "alphaTpad", corners["x1"], corners["y1"], corners["x2"], corners["y2"])
        pad.cd()
        pad.SetTickx()
        pad.SetTicky()

        showAlphaTMet = showAlphaTMet=self.showAlphaTMet and not self.prettyMode

        title = ";"
        if showAlphaTMet :
            title +="#color[%d]{MET/HT}              "%r.kGreen
        title+= "#color[%d]{MHT/HT};#DeltaHT/HT"%r.kBlue
        alphaTHisto = r.TH2D("alphaTHisto",title,100,0.0,1.0,100,0.0,1.0)
        alphaTMetHisto = alphaTHisto.Clone("alphaTMetHisto")

        mht = eventVars["%sSumP4%s"%self.jets].pt() if eventVars["%sSumP4%s"%self.jets] else 0
        met = eventVars[self.met].pt()
        ht  = eventVars["%sSumPt%s"%self.jets]
        deltaHt   = eventVars[self.deltaHtName]
        alphaT    = eventVars["%sAlphaTEt%s"%self.jets]    #hack: hard-coded Et
        alphaTMet = eventVars["%sAlphaTMetEt%s"%self.jets] #hack: hard-coded Et
        
        if ht : alphaTHisto.Fill(mht/ht,deltaHt/ht)
        alphaTHisto.SetStats(False)
        alphaTHisto.SetMarkerStyle(29)
        alphaTHisto.GetYaxis().SetTitleOffset(1.15)
        alphaTHisto.SetMarkerColor(r.kBlue)
        alphaTHisto.GetXaxis().SetTitleSize(self.titleSizeFactor*alphaTHisto.GetXaxis().GetTitleSize())
        alphaTHisto.GetYaxis().SetTitleSize(self.titleSizeFactor*alphaTHisto.GetYaxis().GetTitleSize())
        alphaTHisto.Draw("p")

        if showAlphaTMet :
            if ht : alphaTMetHisto.Fill(met/ht,deltaHt/ht)
            alphaTMetHisto.SetStats(False)
            alphaTMetHisto.SetMarkerStyle(29)
            alphaTMetHisto.GetYaxis().SetTitleOffset(1.25)
            alphaTMetHisto.SetMarkerColor(r.kGreen)
            alphaTMetHisto.Draw("psame")

        legend1 = r.TLegend(0.1, 0.6, 1.0, 0.9)
        legend1.SetBorderSize(0)
        legend1.SetFillStyle(0)
        
        for func in self.alphaFuncs :
            func.Draw("same")
            legend1.AddEntry(func,func.GetName(),"l")

        legend1.Draw()

        #legend2 = r.TLegend(0.1, 0.4, 1.0, 0.6)
        #legend2.SetBorderSize(0)
        #legend2.SetFillStyle(0)
        #legend2.AddEntry(alphaTHisto,"this event","p")
        #if showAlphaTMet :
        #    legend2.AddEntry(alphaTMetHisto,"this event","p")
        #legend2.Draw()
        
        self.canvas.cd()
        pad.Draw()
        return [pad, alphaTHisto, alphaTMetHisto, legend1]

    def rhoPhiPad(self, eventVars, coords, corners):
        pad = r.TPad("rhoPhiPad", "rhoPhiPad", corners["x1"], corners["y1"], corners["x2"], corners["y2"])
        pad.cd()

        skeletonColor = r.kYellow+1
        self.drawSkeleton(coords, skeletonColor)
        self.drawScale(color=skeletonColor, size=0.03, scale=coords["scale"],
                       point={"x":0.0, "y":coords["radius"]+coords["y0"]+0.03})
        return pad

    def drawObjects(self, eventVars=None, etaPhiPad=None, rhoPhiPad=None, rhoPhiCoords=None):
        defArrowSize=0.5*self.arrow.GetDefaultArrowSize()
        defWidth=1

        if self.genMet and not eventVars["isRealData"]:
            rhoPhiPad.cd()
            self.drawGenMet(eventVars, rhoPhiCoords, r.kMagenta, defWidth, defArrowSize*2/6.0)

        arrowSize = defArrowSize
        if not eventVars["isRealData"]:
            for indices, color, printList in self.genParticleIndices:
                rhoPhiPad.cd()
                self.drawGenParticles(eventVars=eventVars,
                                      indices=indices,
                                      coords=rhoPhiCoords,
                                      lineColor=color,
                                      arrowSize=arrowSize)
                arrowSize *= 0.8

                etaPhiPad.cd()
                self.drawGenParticles(eventVars=eventVars,
                                      indices=indices,
                                      lineColor=color,
                                      circleRadius=0.15)

        for jets, indices, color, style, radius in self.jetIndices:
            rhoPhiPad.cd()
            self.drawJets(eventVars=eventVars,
                          jets=jets,
                          indices=indices,
                          coords=rhoPhiCoords,
                          lineColor=color,
                          lineStyle=style,
                          arrowSize=arrowSize,
                          )
            arrowSize *= 0.8
            etaPhiPad.cd()
            self.drawJets(eventVars=eventVars,
                          jets=jets,
                          indices=indices,
                          lineColor=color,
                          lineStyle=style,
                          circleRadius=radius,
                          )

        rhoPhiPad.cd()
        coords = rhoPhiCoords
        if self.jets and self.ra1Mode:
            if not self.prettyMode:
                self.drawHt         (eventVars, coords,r.kBlue+3  , defWidth, defArrowSize*1/6.0)
                self.drawNJetDeltaHt(eventVars, coords,r.kBlue-9  , defWidth, defArrowSize*1/6.0)
            self.drawMht(eventVars, coords,r.kRed     , defWidth, defArrowSize*3/6.0)

        items = [("met", r.kGreen),
                 ("muons", r.kYellow),
                 ("electrons", r.kOrange+7),
                 ("photons", r.kOrange),
                 ] + ([] if self.prettyMode else [("taus", r.kYellow)])

        for item, color in items:
            fixes = getattr(self, item)
            if type(fixes) is str:
                self.legendFunc(lineColor=color, name=fixes, desc=fixes)
                self.drawP4(coords, eventVars[fixes], color, defWidth, arrowSize=defArrowSize*2/6.0)
            elif type(fixes) is tuple and len(fixes) == 2:
                p4Name = "%sP4%s" % fixes
                self.legendFunc(lineColor=color, name=p4Name, desc=p4Name)
                p4Vector = eventVars[p4Name]
                for i in range(len(p4Vector)):
                    self.drawP4(coords, p4Vector.at(i), color, defWidth, arrowSize=defArrowSize*2/6.0)

        if self.recHits and not self.prettyMode:
            #self.drawCleanedRecHits(eventVars, coords,r.kOrange-6, defWidth, defArrowSize*2/6.0)
            self.drawCleanedRecHitSumP4(eventVars, coords,r.kOrange-6, defWidth, defArrowSize*2/6.0)


    def drawLegend(self, corners) :
        pad = r.TPad("legendPad", "legendPad", corners["x1"], corners["y1"], corners["x2"], corners["y2"])
        pad.cd()
        
        legend = r.TLegend(0.0, 0.0, 1.0, 1.0)
        for color, style, desc, gopts in self.legendList:
            self.line.SetLineColor(color)
            self.line.SetLineStyle(style)
            someLine = self.line.DrawLine(0.0, 0.0, 0.0, 0.0)
            legend.AddEntry(someLine, desc, gopts)
        legend.Draw("same")
        self.canvas.cd()
        pad.Draw()
        return [pad,legend]

    def printAllText(self, eventVars, corners) :
        pad = r.TPad("textPad", "textPad", corners["x1"], corners["y1"], corners["x2"], corners["y2"])
        pad.cd()

        defaults = {}
        defaults["size"] = 0.035
        defaults["font"] = 80
        defaults["color"] = r.kBlack
        defaults["slope"] = 0.017
        s = defaults["slope"]

        smaller = {}
        smaller.update(defaults)
        smaller["size"] = 0.034
        
        yy = 0.98
        x0 = 0.01
        x1 = 0.51
        self.printEvent(   eventVars, params = defaults, coords = {"x":x0, "y":yy})

        if self.printExtraText :
            if self.jets:
                self.printJets(eventVars, params = smaller, coords = {"x":x0, "y":yy-2*s}, jets = self.jets, nMax = 7)

            if self.printGen:
                self.printGenJets(  eventVars, params = defaults, coords = {"x":x0, "y":yy-13*s}, nMax = 7)
                self.printGenParticles(eventVars,params=defaults, coords = {"x":x0+0.40, "y":yy-13*s}, nMax = 7)
            if self.jetsOtherAlgo :
                self.printJets(     eventVars, params = smaller, coords = {"x":x0,  "y":yy-13*s}, jets = self.jetsOtherAlgo, nMax = 7)

            if self.muons :
                self.printMuons(    eventVars, params = defaults, coords = {"x":x0, "y":yy-35*s}, muons = self.muons, nMax = 3)
            self.printVertices(     eventVars, params = defaults, coords = {"x":x1, "y":yy-35*s}, nMax = 3)

            if self.photons :
                self.printPhotons(  eventVars, params = defaults, coords = {"x":x0, "y":yy-42*s}, photons = self.photons, nMax = 3)
            if self.electrons :
                self.printElectrons(eventVars, params = defaults, coords = {"x":x1, "y":yy-42*s}, electrons = self.electrons, nMax = 3)

            if not self.prettyMode :
                if self.recHits :
                    self.printRecHits(eventVars,params = defaults,coords = {"x":x0, "y":yy-49*s}, recHits = self.recHits, nMax = 3)
                if self.recHitsOtherAlgo :
                    self.printRecHits(eventVars,params = defaults,coords = {"x":x1, "y":yy-49*s}, recHits = self.recHitsOtherAlgo, nMax = 3)

            #if self.flagsToPrint :
            #    self.printFlags(    eventVars, params = defaults, coords = {"x":x0, "y":yy-49*s}, flags = self.flagsToPrint)
            if self.ra1Mode :
                self.printKinematicVariables(eventVars, params = defaults, coords = {"x":x0, "y":yy-25*s}, jets = self.jets, jets2 = self.jetsOtherAlgo)
                if self.ra1CutBits :
                    self.printCutBits(       eventVars, params = defaults, coords = {"x":x0, "y":yy-30*s}, jets = self.jets, jets2 = self.jetsOtherAlgo,
                                         met = self.met, met2 = self.metOtherAlgo)
        self.canvas.cd()
        pad.Draw()
        return [pad]

    def display(self, eventVars) :
        rhoPhiPadYSize = 0.50*self.canvas.GetAspectRatio()
        rhoPhiPadXSize = 0.50
        radius = 0.4

        rhoPhiCoords = {"scale":self.scale, "radius":radius,
                        "x0":radius, "y0":radius+0.05}

        rhoPhiCorners = {"x1":0.0,
                         "y1":0.0,
                         "x2":rhoPhiPadXSize,
                         "y2":rhoPhiPadYSize}

        etaPhiCorners = {"x1":rhoPhiPadXSize - 0.18,
                         "y1":rhoPhiPadYSize - 0.08*self.canvas.GetAspectRatio(),
                         "x2":rhoPhiPadXSize + 0.12,
                         "y2":rhoPhiPadYSize + 0.22*self.canvas.GetAspectRatio()}

        alphaCorners = {"x1":rhoPhiPadXSize - 0.08,
                        "y1":0.0,
                        "x2":rhoPhiPadXSize + 0.12,
                        "y2":0.55}

        legendCorners = {"x1":0.0,
                         "y1":rhoPhiPadYSize,
                         "x2":1.0-rhoPhiPadYSize,
                         "y2":1.0}

        textCorners =  {"x1":rhoPhiPadXSize + 0.11,
                        "y1":0.0,
                        "x2":1.0,
                        "y2":1.0}

        rhoPhiPad = self.rhoPhiPad(eventVars, rhoPhiCoords, rhoPhiCorners)
        etaPhiPad, etaPhiPlot = self.etaPhiPad(eventVars, etaPhiCorners)

        keep = [rhoPhiPad, etaPhiPad, etaPhiPlot]
        self.drawObjects(eventVars, etaPhiPad, rhoPhiPad, rhoPhiCoords)

        self.canvas.cd()
        rhoPhiPad.Draw()

        if self.ra1Mode:
            keep.append(self.ra1EtaPhi(eventVars, etaPhiPad))
            keep.append(self.drawAlphaPlot(eventVars, r.kBlack, alphaCorners))

        etaPhiPad.Draw()

        keep.append(self.drawLegend(corners=legendCorners))
        keep.append(self.printAllText(eventVars, corners=textCorners))
        return keep
