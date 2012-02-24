import os,collections,copy,ROOT as r
from supy import analysisStep,utils,configuration
#####################################
pdgLookupExists = False
try:
    import pdgLookup
    pdgLookupExists = True
except ImportError:
    pass
#####################################
class ttbar(analysisStep) :
    
    def __init__(self, jets = None, met = None, muons = None, electrons = None, photons = None, taus = None, scale = 200 , prettyMode = True , printExtraText = False ) :

        for item in ["jets","met","muons","electrons","photons","taus","scale","prettyMode","printExtraText"] :
            setattr(self,item,eval(item))

        self.etaBE = configuration.detectorSpecs()["cms"]["etaBE"]
        self.jetRadius = 0.7 if "7" in self.jets[0] else 0.5
        self.genJets = "gen%sGenJetsP4"%(self.jets[0].replace("xc","")[:3])
        self.genMet  = "genmetP4True"
        self.titleSizeFactor = 1.0
        self.legendDict = collections.defaultdict(int)
        self.legendList = []
        
        self.prettyReName = {
            "clean jets (xcak5JetPFPat)": "jets (AK5 PF)",
            "MET (metP4PF)": "PF MET",
            "muons (muonPF)": "muons",
            "electrons (electronPF)": "electrons",
            "xcak5JetPFPat": "AK5 PF Jets",
            "muonPF": "muons",
            "electronPF": "electrons",
            }

    def outputSuffix(self) :
        return "_displaysTTbar.root"
    
    def uponAcceptance(self, eV) :
        self.keep = []

        r.gStyle.SetOptStat(110011)

        self.preparePads(eV)
        self.drawMet(eV, r.kRed, 1, self.arrow.GetDefaultArrowSize() )
        self.drawGenTTdecay(eV)

        self.drawRhoPhiPlot(eV)
        self.drawEtaPhiPlot(eV)
        self.drawLegend()
        self.printNarrowText(eV)
        
        self.outputFile.cd()
        self.canvas.Write("canvas_%d"%self.canvasIndex)
        self.canvasIndex += 1
        del self.keep

    def setup(self, chain, fileDir) :
        someDir = r.gDirectory
        self.outputFile = r.TFile(self.outputFileName, "RECREATE")
        someDir.cd()
        
        self.ellipse = r.TEllipse(); self.ellipse.SetFillStyle(0)
        self.line = r.TLine()
        self.arrow = r.TArrow(); self.arrow.SetDefaultArrowSize( 0.6 * self.arrow.GetDefaultArrowSize() )
        self.text = r.TText()
        self.latex = r.TLatex()
        self.canvas = utils.canvas("canvas"); self.canvas.SetFixedAspectRatio(); self.canvasIndex = 0

        rhoPhiPadYSize = 0.40*self.canvas.GetAspectRatio()
        rhoPhiPadXSize = 0.40
        radius = 0.5
        offset = 0.02

        self.rhoPhiCoords = {"scale":self.scale, "radius":radius, "x0":radius, "y0":radius}

        rhoPhiCorners = {"x1":offset,                   "x2":offset + rhoPhiPadXSize,    "y1":1-offset - rhoPhiPadYSize, "y2":1-offset}
        etaPhiCorners = {"x1":offset + rhoPhiPadXSize,  "x2":offset + 2*rhoPhiPadXSize,  "y1":1-offset - rhoPhiPadYSize, "y2":1-offset}
        legendCorners = {"x1":offset,                   "x2":offset + 0.5*rhoPhiPadXSize,"y1":offset,                    "y2":1-offset*2 - rhoPhiPadYSize}
        narrowCorners = {"x1": 2*rhoPhiPadXSize,        "x2":1.0,                        "y1":offset,                    "y2":1.0 - offset}

        self.etaPhiPad = r.TPad("etaPhiPad", "etaPhiPad", etaPhiCorners['x1'], etaPhiCorners['y1'], etaPhiCorners['x2'], etaPhiCorners['y2'] )
        self.rhoPhiPad = r.TPad("rhoPhiPad", "rhoPhiPad", rhoPhiCorners['x1'], rhoPhiCorners['y1'], rhoPhiCorners['x2'], rhoPhiCorners['y2'] )
        self.legendPad = r.TPad("legendPad", "legendPad", legendCorners["x1"], legendCorners["y1"], legendCorners["x2"], legendCorners["y2"] )
        self.narrowPad = r.TPad("narrowPad", "narrowPad", narrowCorners["x1"], narrowCorners["y1"], narrowCorners["x2"], narrowCorners["y2"] )

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

    def printEvent(self, eV, params, coords) :
        self.prepareText(params, coords)
        for message in ["Run   %#10d"%eV["run"],
                        "Ls    %#10d"%eV["lumiSection"],
                        "Event %#10d"%eV["event"],
                        "genQ (GeV) %#5.1f"%eV["genQ"] if not eV["isRealData"] else "",
                        "ttDecayMode:  %s"%eV["ttDecayChannel"],
                        ] :
            if message : self.printText(message)
        
    
    def drawSkeleton(self, color) :
        coords = self.rhoPhiCoords
        r.gPad.AbsCoordinates(False)
        
        self.ellipse.SetLineColor(color)
        self.ellipse.SetLineWidth(1)
        self.ellipse.SetLineStyle(1)
        self.ellipse.DrawEllipse(coords["x0"], coords["y0"], coords["radius"], coords["radius"], 0.0, 360.0, 0.0, "")
        self.ellipse.DrawEllipse(coords["x0"], coords["y0"], coords["radius"]/10, coords["radius"]/10, 0.0, 360.0, 0.0, "")

        self.line.SetLineColor(color)
        self.line.DrawLine(coords["x0"]-coords["radius"], coords["y0"]                 , coords["x0"]+coords["radius"], coords["y0"]                 )
        self.line.DrawLine(coords["x0"]                 , coords["y0"]-coords["radius"], coords["x0"]                 , coords["y0"]+coords["radius"])

        self.drawScale(color = color, size = 0.03, scale = coords["scale"], point = {"x":0.03, "y":coords["radius"]+coords["y0"]+0.03})

    def drawScale(self, color, size, scale, point) :
        self.latex.SetTextSize(size)
        self.latex.SetTextColor(color)
        self.latex.DrawLatex(point["x"], point["y"],"radius = "+str(scale)+" GeV p_{T}")

    def drawP4(self, c, p4, color, lineWidth, arrowSize, p4Initial = None) :
        x0 = c["x0"]+p4Initial.px()*c["radius"]/c["scale"] if p4Initial else c["x0"]
        y0 = c["y0"]+p4Initial.py()*c["radius"]/c["scale"] if p4Initial else c["y0"]
        x1 = x0+p4.px()*c["radius"]/c["scale"]
        y1 = y0+p4.py()*c["radius"]/c["scale"]

        self.arrow.SetLineColor(color)
        self.arrow.SetLineWidth(lineWidth)
        self.arrow.SetArrowSize(arrowSize*p4.pt()/c["scale"]*2)
        self.arrow.SetFillColor(color)
        self.arrow.DrawArrow(x0,y0,x1,y1)
        
    def drawCircle(self, p4, color, lineWidth, circleRadius, lineStyle = 1) :
        self.ellipse.SetLineColor(color)
        self.ellipse.SetLineWidth(lineWidth)
        self.ellipse.SetLineStyle(lineStyle)
        self.ellipse.DrawEllipse(p4.eta(), p4.phi(), circleRadius, circleRadius, 0.0, 360.0, 0.0, "")
            
    def drawGenParticles(self, eV, color, arrow, circle,  status = [], pdgs = [], moms = [], label = "") :
        self.legendFunc(color, name = "genParticle"+label, desc = label)
        p4 = eV["genP4"]
        for iParticle in range(6,len(p4)) :
            particle = p4[iParticle]
            if ( eV["genStatus"].at(iParticle)   not in status or
                 eV["genPdgId"].at(iParticle)    not in pdgs or
                 eV["genMotherPdgId"][iParticle] not in moms) : continue

            self.rhoPhiPad.cd(); self.drawP4(self.rhoPhiCoords, particle, color, arrow['width'], arrow['size'])
            self.etaPhiPad.cd(); self.drawCircle( particle, color, circle['width'], circle['radius'])
            
    def drawGenTTdecay(self, eV) :
        if not eV["isRealData"] :
            arrow = {"width":1, "size":0.5*self.arrow.GetDefaultArrowSize()}
            circle = {"width":3, "radius":0.08}
            status = [3]
            self.drawGenParticles( eV,            9, arrow, circle, status, pdgs = [5,-5],                 moms = [6,-6],           label = "b from top")
            self.drawGenParticles( eV,            7, arrow, circle, status, pdgs = [-4,-3,-2,-1,1,2,3,4],  moms = [24,-24],         label = "quark from W")
            self.drawGenParticles( eV,    r.kSpring, arrow, circle, status, pdgs = [-15,-13,-11,11,13,15], moms = [24,-24],         label = "lepton from W")
            self.drawGenParticles( eV, r.kOrange+1,  arrow, circle, status, pdgs = [-16,-14,-12,12,14,16], moms = [24,-24],         label = "neutrino from W")
            self.drawGenParticles( eV,           28, arrow, circle, status, pdgs = [21],                   moms = range(-6,7)+[21], label = "gluon (status 3)")

    def drawMet(self, eV, color, lineWidth, arrowSize) :
        if not self.met: return
        self.legendFunc(color, name = "met%s"%self.met, desc = "MET (%s)"%self.met)
        self.line.SetLineColor(color)
        self.rhoPhiPad.cd();  self.drawP4(self.rhoPhiCoords, eV[self.met], color, lineWidth, arrowSize)
        self.etaPhiPad.cd();  self.line.DrawLine( -3, eV[self.met].phi(), 3, eV[self.met].phi(),  )

    def drawCleanJets(self, eV, coords, jets, color, lineWidth, arrowSize) :
        self.legendFunc(color, name = "cleanJet".join(jets), desc = "clean jets (%s%s)"%jets)
        
        p4s = eV['CorrectedP4'.join(jets)]
        
        cleanJetIndices = eV["Indices".join(jets)]
        for iJet in cleanJetIndices :
            self.drawP4(coords, p4s.at(iJet), color, lineWidth, arrowSize)
                    
    def drawMuons(self, eV, coords, color, lineWidth, arrowSize) :
        self.legendFunc(color, name = "%smuon%s"%self.muons, desc = "muons (%s%s)"%self.muons)
        p4Vector=eV["%sP4%s"%self.muons]
        for i in range(len(p4Vector)) :
            self.drawP4(coords, p4Vector.at(i), color, lineWidth, arrowSize)
            
    def drawElectrons(self, eV, coords, color, lineWidth, arrowSize) :
        self.legendFunc(color, name = "%selectron%s"%self.electrons, desc = "electrons (%s%s)"%self.electrons)
        p4Vector=eV["%sP4%s"%self.electrons]
        for i in range(len(p4Vector)) :
            self.drawP4(coords, p4Vector.at(i), color, lineWidth, arrowSize)
            
    def preparePads(self, eV) :
        self.legendPad.Clear();
        self.rhoPhiPad.Clear();
        self.etaPhiPad.Clear();
        
        self.rhoPhiPad.cd()
        self.drawSkeleton(r.kYellow+1)

        self.etaPhiPad.cd(); 
        self.etaPhiPad.SetTickx(); self.etaPhiPad.SetTicky()
        etaPhiPlot = r.TH2D("etaPhi",";#eta;#phi;",1, -3.0, 3.0, 1, -r.TMath.Pi(), r.TMath.Pi() )
        etaPhiPlot.SetStats(False); etaPhiPlot.SetTitle(""); etaPhiPlot.Draw()
        self.keep.append(etaPhiPlot)

        self.line.SetLineColor(r.kBlack)
        self.line.DrawLine(-self.etaBE, etaPhiPlot.GetYaxis().GetXmin(), -self.etaBE, etaPhiPlot.GetYaxis().GetXmax() )
        self.line.DrawLine( self.etaBE, etaPhiPlot.GetYaxis().GetXmin(),  self.etaBE, etaPhiPlot.GetYaxis().GetXmax() )

    def drawEtaPhiPlot (self, eV) :
        self.etaPhiPad.cd()
        
        jets = eV["CorrectedP4".join(self.jets)]
        for index in eV["Indices".join(self.jets)] :
            jet = jets.at(index)
            self.drawCircle(jet, r.kBlue, lineWidth = 1, circleRadius = self.jetRadius)
            if jet.pt()>35 :
                self.drawCircle(jet, r.kBlue, lineWidth = 1, circleRadius = self.jetRadius*(1 + 0.0005*(jet.pt()-30)), lineStyle=1)
                self.drawCircle(jet, r.kBlue, lineWidth = 1, circleRadius = self.jetRadius*(1 - 0.0005*(jet.pt()-30)), lineStyle=1)
        
        self.canvas.cd()
        self.etaPhiPad.Draw()

    def drawRhoPhiPlot(self, eV):
        coords = self.rhoPhiCoords

        self.rhoPhiPad.cd()
        defArrowSize = self.arrow.GetDefaultArrowSize()
        defWidth=1

        if self.jets :      self.drawCleanJets  (eV, coords, self.jets, r.kBlue  , defWidth, defArrowSize)
        if self.muons :     self.drawMuons      (eV, coords,            r.kGreen+3 , defWidth, defArrowSize)
        if self.electrons : self.drawElectrons  (eV, coords,            r.kGreen+3, defWidth, defArrowSize)

        self.canvas.cd()
        self.rhoPhiPad.Draw()

    def drawLegend(self) :
        self.legendPad.cd();
        
        legend = r.TLegend(0.0, 0.0, 1.0, 1.0)
        for item in self.legendList :
            self.line.SetLineColor(item[0])
            someLine = self.line.DrawLine(0.0,0.0,0.0,0.0)
            legend.AddEntry(someLine, item[1], item[2])
        legend.Draw("same")
        self.canvas.cd()
        self.legendPad.Draw()
        self.keep.append(legend)

    def renamedDesc(self, desc) :
        if not self.prettyMode : return desc
        elif desc in self.prettyReName : return self.prettyReName[desc]
        else : return desc
        
    def legendFunc(self, color, name, desc) :
        if not self.legendDict[name] :
            self.legendDict[name] = True
            self.legendList.append( (color, self.renamedDesc(desc), "l") )

    def printAllText(self, eV, corners) :
        pad = r.TPad("textPad", "textPad", corners["x1"], corners["y1"], corners["x2"], corners["y2"])
        pad.cd()

        defaults = {}
        defaults["size"] = 0.035
        defaults["font"] = 80
        defaults["color"] = r.kBlack
        defaults["slope"] = 0.017
        s = defaults["slope"]

        smaller = copy.copy(defaults)
        smaller["size"] = 0.034
        
        yy = 0.98
        x0 = 0.01
        x1 = 0.45

        if self.printExtraText :
            self.printVertices(eV, params = defaults, coords = {"x":x1, "y":yy}, nMax = 3)
            self.printJets(    eV, params = defaults, coords = {"x":x0, "y":yy-7*s}, jets = self.jets, nMax = 7)

            if self.doGenJets :
                self.printGenJets(  eV, params = defaults, coords = {"x":x0,      "y":yy-18*s}, nMax = 7)
                self.printGenParticles(eV,params=defaults, coords = {"x":x0+0.40, "y":yy-18*s}, nMax = 7)
            if self.photons :
                self.printPhotons(  eV, params = defaults, coords = {"x":x0,      "y":yy-40*s}, photons = self.photons, nMax = 3)
            if self.electrons :
                self.printElectrons(eV, params = defaults, coords = {"x":x0+0.50, "y":yy-40*s}, electrons = self.electrons, nMax = 3)
            if self.muons :
                muonPars = defaults if self.prettyMode else smaller
                self.printMuons(    eV, params = muonPars, coords = {"x":x0,      "y":yy-47*s}, muons = self.muons, nMax = 3)

        self.canvas.cd()
        pad.Draw()
        return [pad]

    def printNarrowText(self, eV) :
        self.narrowPad.Clear(); self.narrowPad.cd()
        defaults = {"size":0.085, "font":80, "color":r.kBlack, "slope":0.026}
        self.printEvent(   eV, params = defaults, coords = {"x":0.01, "y":0.98})
        self.canvas.cd()
        self.narrowPad.Draw()

    def printVertices(self, eV, params, coords, nMax) :
        self.prepareText(params, coords)
        self.printText("Vertices")
        self.printText("ID   Z(cm)%s"%(" sumPt(GeV)" if not self.prettyMode else ""))
        self.printText("----------%s"%("-----------" if not self.prettyMode else ""))

        nVertices = eV["vertexNdof"].size()
        for i in range(nVertices) :
            if nMax<=i :
                self.printText("[%d more not listed]"%(nVertices-nMax))
                break
            
            out = "%2s  %6.2f"%("G " if i in eV["vertexIndices"] else "  ", eV["vertexPosition"].at(i).z())
            if not self.prettyMode : out += " %5.0f"%eV["vertexSumPt"].at(i)
            self.printText(out)

    def printElectrons(self, eV, params, coords, electrons, nMax) :
        self.prepareText(params, coords)
        p4Vector = eV["%sP4%s"        %electrons]
        cIso = eV["%sIsoCombined%s"%electrons]
        ninetyFive = eV["%sID95%s"%electrons]
     
        self.printText(self.renamedDesc(electrons[0]+electrons[1]))
        self.printText("ID   pT  eta  phi  cIso")
        self.printText("-----------------------")

        nElectrons = p4Vector.size()
        for iElectron in range(nElectrons) :
            if nMax<=iElectron :
                self.printText("[%d more not listed]"%(nElectrons-nMax))
                break
            electron=p4Vector[iElectron]

            outString = "%2s"%("95" if ninetyFive[iElectron] else "  ")
            outString+="%5.0f %4.1f %4.1f"%(electron.pt(), electron.eta(), electron.phi())
            outString+=" %5.2f"%cIso[iElectron] if cIso[iElectron]!=None else " %5s"%"-"
            self.printText(outString)

    def printMuons(self, eV, params, coords, muons, nMax) :
        self.prepareText(params, coords)
        p4Vector = eV["%sP4%s"     %muons]
        tight    = eV["%sIDtight%s"%muons]
        iso      = eV["%sCombinedRelativeIso%s"%muons]
        tr       = eV["%sIsTrackerMuon%s"%muons]
        gl       = eV["%sIsGlobalMuon%s"%muons]
        glpt     = eV["%sIDGlobalMuonPromptTight%s"%muons]
        
        self.printText(self.renamedDesc(muons[0]+muons[1]))
        self.printText("ID   pT  eta  phi  cIso cat")
        self.printText("---------------------------")

        nMuons = p4Vector.size()
        for iMuon in range(nMuons) :
            if nMax<=iMuon :
                self.printText("[%d more not listed]"%(nMuons-nMax))
                break
            muon=p4Vector[iMuon]

            outString = "%1s%1s"% (" ","T" if tight[iMuon] else " ")
            outString+= "%5.0f %4.1f %4.1f"%(muon.pt(), muon.eta(), muon.phi())
            outString+= " %5.2f"%(iso[iMuon]) if iso[iMuon]<100.0 else ">100".rjust(6)
            outString+= " %s%s%s"%("T" if tr[iMuon] else " ", "G" if gl[iMuon] else " ","P" if glpt[iMuon] else " ")

            self.printText(outString)

    def printJets(self, eV, params, coords, jets, nMax) :
        self.prepareText(params, coords)
        jets2 = (jets[0].replace("xc",""),jets[1])
        isPf = "PF" in jets[0]
        
        p4Vector         = eV['%sCorrectedP4%s'%jets]
        corrVector       = eV['%sCorrFactor%s'      %jets2]

        if not isPf :
            jetEmfVector  = eV['%sEmEnergyFraction%s'%jets2]
            jetFHpdVector = eV['%sJetIDFHPD%s'       %jets2]
            jetFRbxVector = eV['%sJetIDFRBX%s'       %jets2]
            jetN90Vector  = eV['%sJetIDN90Hits%s'    %jets2]
            
            loose = eV["%sJetIDloose%s"%jets2]
            tight = eV["%sJetIDtight%s"%jets2]
            
        else :
            chf = eV["%sFchargedHad%s"%jets2]
            nhf = eV["%sFneutralHad%s"%jets2]

            cef = eV["%sFchargedEm%s"%jets2]
            nef = eV["%sFneutralEm%s"%jets2]

            cm  = eV["%sNcharged%s"%jets2]
            nm  = eV["%sNneutral%s"%jets2]
            
            loose = eV["%sPFJetIDloose%s"%jets2]
            tight = eV["%sPFJetIDtight%s"%jets2]
            
        self.printText(self.renamedDesc(jets[0]+jets[1]))
        self.printText("ID   pT  eta  phi%s"%("   EMF  fHPD  fRBX N90 corr" if not isPf else "   CHF  NHF  CEF  NEF CM corr"))
        self.printText("-----------------%s"%("---------------------------" if not isPf else "-----------------------------"))

        nJets = p4Vector.size()
        for iJet in range(nJets) :
            if nMax<=iJet :
                self.printText("[%d more not listed]"%(nJets-nMax))
                break
            jet=p4Vector[iJet]

            outString = "%1s%1s"% ("L" if loose[iJet] else " ", "T" if tight[iJet] else " ")
            outString+="%5.0f %4.1f %4.1f"%(jet.pt(), jet.eta(), jet.phi())

            if not isPf :
                outString+=" %5.2f %5.2f %5.2f %3d %4.2f"%(jetEmfVector.at(iJet), jetFHpdVector.at(iJet), jetFRbxVector.at(iJet), jetN90Vector.at(iJet), corrVector.at(iJet))
            else :
                outString+=" %5.3f %4.2f %4.2f %4.2f%3d %4.2f"%(chf.at(iJet), nhf.at(iJet), cef.at(iJet), nef.at(iJet), cm.at(iJet), corrVector.at(iJet))
            self.printText(outString)

    def printGenJets(self, eV, params, coords, nMax) :
        self.prepareText(params, coords)
        p4Vector = eV[self.genJets]
            
        self.printText(self.renamedDesc(self.genJets))
        self.printText("   pT  eta  phi")
        self.printText("---------------")

        nJets = p4Vector.size()
        for iJet in range(nJets) :
            if nMax<=iJet :
                self.printText("[%d more not listed]"%(nJets-nMax))
                break
            jet = p4Vector[iJet]
            self.printText("%5.0f %4.1f %4.1f"%(jet.pt(), jet.eta(), jet.phi()))

    def printGenParticles(self, eV, params, coords, nMax) :
        self.prepareText(params, coords)
        p4s    = eV["genP4"]
        status = eV["genStatus"]
        ids    = eV["genPdgId"]
        
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


    def mergeFunc(self, products) :
        def psFromRoot(listOfInFileNames, outFileName) :
            if not len(listOfInFileNames) : return
            options = "pdf"
            dummyCanvas = utils.canvas("display")
            dummyCanvas.Print(outFileName+"[", options)
            for inFileName in listOfInFileNames :
                inFile = r.TFile(inFileName)
                keys = inFile.GetListOfKeys()
                for key in keys :
                    someObject = inFile.Get(key.GetName())
                    if someObject.ClassName()!="TCanvas" : print "Warning: found an object which is not a TCanvas in the display root file"
                    someObject.Print(outFileName, options)
                inFile.Close()
                os.remove(inFileName)                    
            dummyCanvas.Print(outFileName+"]", options)
            print "The display file \""+outFileName+"\" has been written."    
        
        psFromRoot(products["outputFileName"], self.outputFileName.replace(".root", ".pdf"))
        print utils.hyphens

    def endFunc(self, chains) :
        self.outputFile.Write(); self.outputFile.Close(); del self.canvas

