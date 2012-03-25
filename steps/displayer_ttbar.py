import math,os,collections,configuration, ROOT as r,numpy as np
import supy
#####################################
class ttbar(supy.steps.displayer) :
    
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
        self.keep = []
        
        self.prettyReName = {
            "clean jets (xcak5JetPFPat)": "jets (AK5 PF)",
            "MET (metP4PF)": "PF MET",
            "muons (muonPF)": "muons",
            "electrons (electronPF)": "electrons",
            "xcak5JetPFPat": "AK5 PF Jets",
            "muonPF": "muons",
            "electronPF": "electrons",
            }

    def reset(self) :
        del self.keep
        self.keep = []
        self.preparePads()

    def display(self, eV) :
        r.gStyle.SetOptStat(110011)

        self.drawMet(eV, r.kRed+1, 1 )
        self.drawJets(eV, r.kBlue+3, 1 )
        self.drawLeptons(eV, r.kGreen+3, 1, kind = "muon" , desc = 'charged lepton')
        self.drawLeptons(eV, r.kGreen+3, 1, kind = "electron", desc = 'charged lepton' )
        self.drawTopReco(eV)
        self.drawGenTTdecay(eV)

        self.canvas.cd()
        self.etaPhiPad.Draw()
        self.rhoPhiPad.Draw()

        self.drawLegend()
        self.printNarrowText(eV)

    def setup(self, chain, fileDir) :
        super(ttbar, self).setup(chain,fileDir)
        
        self.ellipse = r.TEllipse(); self.ellipse.SetFillStyle(0)
        self.metunc= r.TEllipse(); self.metunc.SetFillStyle(0); self.metunc.SetLineStyle(3); self.metunc.SetLineColor(r.kRed+1)
        self.line = r.TLine()
        self.arrow = r.TArrow(); self.arrow.SetDefaultArrowSize( 0.6 * self.arrow.GetDefaultArrowSize() )
        self.marker = r.TMarker();
        self.text = r.TText()
        self.latex = r.TLatex()

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

    def scaleRho(self,gev,coords=None) :
        if coords==None : coords = self.rhoPhiCoords
        return gev * coords['radius']/coords['scale']

    def drawP4(self, c, p4, color, lineWidth, arrowSize, p4Initial = None) :
        self.rhoPhiPad.cd()
        x0 = c["x0"]+self.scaleRho(p4Initial.px()) if p4Initial else c["x0"]
        y0 = c["y0"]+self.scaleRho(p4Initial.py()) if p4Initial else c["y0"]
        x1 = x0+self.scaleRho(p4.px())
        y1 = y0+self.scaleRho(p4.py())

        self.arrow.SetLineColor(color)
        self.arrow.SetLineWidth(lineWidth)
        self.arrow.SetArrowSize(arrowSize*p4.pt()/c["scale"]*2)
        self.arrow.SetFillColor(color)
        self.arrow.DrawArrow(x0,y0,x1,y1)
        
    def drawCircle(self, p4, color, lineWidth, circleRadius, lineStyle = 1) :
        self.etaPhiPad.cd()
        self.ellipse.SetLineColor(color)
        self.ellipse.SetLineWidth(lineWidth)
        self.ellipse.SetLineStyle(lineStyle)
        self.ellipse.DrawEllipse(p4.eta(), p4.phi(), circleRadius, circleRadius, 0.0, 360.0, 0.0, "")

    def drawMarker(self, p4, color, size, style, ptMode=False) :
        self.marker.SetMarkerColor(color)
        self.marker.SetMarkerSize(size)
        self.marker.SetMarkerStyle(style)
        phi = p4.phi()
        if not ptMode :
            self.etaPhiPad.cd()
            self.marker.DrawMarker(p4.eta(), phi)
        else :
            self.rhoPhiPad.cd()
            N = 1 if p4.pt()<200 else 200/p4.pt()
            self.marker.DrawMarker(self.rhoPhiCoords['x0'] + N*self.scaleRho(p4.x()), self.rhoPhiCoords['y0'] + N*self.scaleRho(p4.y()))
            
    def drawGenParticles(self, eV, color, arrow, circle,  status = [], pdgs = [], moms = [], label = "") :
        self.legendFunc(color, name = "genParticle"+label, desc = label)
        p4 = eV["genP4"]
        for iParticle in range(6,len(p4)) :
            particle = p4[iParticle]
            if ( eV["genStatus"].at(iParticle)   not in status or
                 eV["genPdgId"].at(iParticle)    not in pdgs or
                 eV["genMotherPdgId"][iParticle] not in moms) : continue

            self.drawP4(self.rhoPhiCoords, particle, color, arrow['width'], arrow['size'])
            self.drawMarker( particle, color, 0.5, r.kFullCircle if eV["genPdgId"].at(iParticle)>0 else r.kOpenCircle )
            
    def drawGenTTdecay(self, eV) :
        if not eV["isRealData"] :
            lSign = -1 if any( eV["genPdgId"][i] in [11,13,15] for i in range(len(eV["genPdgId"])) if eV["genMotherPdgId"][i]==-24 ) else 1
            arrow = {"width":1, "size":0.5*self.arrow.GetDefaultArrowSize()}
            circle = {"width":3, "radius":0.08}
            status = [3]
            
            self.drawGenParticles( eV, r.kMagenta-6, arrow, circle, status, pdgs = [+lSign*5],              moms = [6,-6],         label = "gen b from top (blv)")
            self.drawGenParticles( eV,    r.kBlue-6, arrow, circle, status, pdgs = [-lSign*5],              moms = [6,-6],         label = "gen b from top (bpq)")
            self.drawGenParticles( eV,            7, arrow, circle, status, pdgs = [-4,-3,-2,-1,1,2,3,4],  moms = [24,-24],        label = "gen quark from W")
            self.drawGenParticles( eV,    r.kSpring, arrow, circle, status, pdgs = [-15,-13,-11,11,13,15], moms = [24,-24],        label = "gen lepton from W")
            self.drawGenParticles( eV,    r.kOrange, arrow, circle, status, pdgs = [-16,-14,-12,12,14,16], moms = [24,-24],        label = "gen neutrino from W")
            self.drawGenParticles( eV,           28, arrow, circle, status, pdgs = [21],                   moms = range(-6,7)+[21],label = "gen gluon")
            p4 = eV['genP4']
            pdg = eV['genPdgId']
            top = next((p4[i] for i in range(len(p4)) if pdg[i] == +6 ), None)
            bar = next((p4[i] for i in range(len(p4)) if pdg[i] == -6 ), None)
            if top and bar :
                self.drawMarker( top, r.kBlack, 1, r.kFullStar)
                self.drawMarker( top, r.kBlack, 1, r.kFullStar, ptMode = True)
                self.drawMarker( bar, r.kBlack, 1, r.kOpenStar)
                self.drawMarker( bar, r.kBlack, 1, r.kOpenStar, ptMode = True)
                self.etaPhiPad.cd()
                self.line.SetLineWidth(1)
                self.line.SetLineColor(r.kBlack)
                self.line.DrawLine( 0, -3.5, top.eta()-bar.eta(), -3.5  )


            scale = 0.001
            self.line.SetLineWidth(5); self.line.SetLineColor(r.kBlack);  self.line.DrawLine( p4[4].pz() * scale, 3.8, p4[5].pz() * scale, 3.8 )
            self.line.SetLineWidth(2); self.line.SetLineColor(r.kWhite);  self.line.DrawLine(             0, 3.8,      p4[4 if pdg[4]<0 else 5].pz() * scale, 3.8 )
            self.line.SetLineWidth(1); self.line.SetLineColor(r.kBlack);  self.line.DrawLine( 0, 3.8, 0, 3.5)


    def drawTopReco(self, eV) :
        def draw(p4, color) :
            self.rhoPhiPad.cd(); self.drawP4(self.rhoPhiCoords, p4, color, 1, self.arrow.GetDefaultArrowSize()*0.6)
            self.etaPhiPad.cd(); self.drawCircle(p4, color, 1, 0.2)
        reco = eV["TopReconstruction"][0]
        draw(reco['lep'], r.kGreen+2)
        draw(reco['nu'], r.kOrange+2)
        draw(reco['hadP'], r.kCyan+2)
        draw(reco['hadQ'], r.kCyan+2)
        draw(reco['hadB'], r.kBlue+2)
        draw(reco['lepB'], r.kMagenta+2)
        if reco['iX']!=None : self.etaPhiPad.cd(); self.drawCircle(eV["CorrectedP4".join(self.jets)][reco['iX']], 48, 1, 0.2)
        self.etaPhiPad.cd()
        self.drawMarker(reco['top'], 40, 1.5, r.kFullStar)
        self.drawMarker(reco['top'], 40, 1.5, r.kFullStar, ptMode = True)
        self.drawMarker(reco['tbar'], 40, 1.5, r.kOpenStar)
        self.drawMarker(reco['tbar'], 40, 1.5, r.kOpenStar, ptMode = True)
        self.etaPhiPad.cd()
        self.line.SetLineWidth(1)
        self.line.SetLineColor(40)
        self.line.DrawLine( 0, -3.6, reco['top'].eta() - reco['tbar'].eta(), -3.6 )
        xp,xm = eV["fitTopPartonXplusminus"]
        scale = 3.5
        self.line.SetLineWidth(2); self.line.SetLineColor(r.kGray);  self.line.DrawLine( scale * xm, 3.6, scale * xp, 3.6 )

        if reco['lepBound'] : self.drawLepBoundCurve(reco['lep'])

    def drawLepBoundCurve(self, lep) :
        self.rhoPhiPad.cd()
        self.marker.SetMarkerColor(r.kGreen)
        self.marker.SetMarkerSize(0.2)
        self.marker.SetMarkerStyle(6)
        lPt = lep.pt()
        lPhi = lep.phi()
        def draw(pt,phi) :
            self.marker.DrawMarker(self.rhoPhiCoords['x0']+self.scaleRho(pt) * math.cos(phi),
                                   self.rhoPhiCoords['y0']+self.scaleRho(pt) * math.sin(phi))
        for i in range(100) :
            dphi = math.pi - 0.1*i
            nuPt = 80.4**2 / (2*lPt*(1-math.cos(dphi)))
            draw( nuPt, lPhi+dphi)
            draw( nuPt, lPhi-dphi)
            

    def drawMet(self, eV, color, lineWidth) :
        if not self.met: return
        self.legendFunc(color, name = "met%s"%self.met, desc = "MET (%s)"%self.met)
        self.line.SetLineColor(color)
        self.etaPhiPad.cd();  self.line.DrawLine( -3, eV[self.met].phi(), 3, eV[self.met].phi()  )
        self.rhoPhiPad.cd();  self.drawP4(self.rhoPhiCoords, eV[self.met], color, lineWidth, self.arrow.GetDefaultArrowSize() )

        coords=self.rhoPhiCoords
        x0 = coords['x0']+self.scaleRho(eV[self.met].px())
        y0 = coords['x0']+self.scaleRho(eV[self.met].py())
        eig,Rinv = np.linalg.eig(eV["metCovariancePF"])
        self.metunc.DrawEllipse(x0,y0,self.scaleRho(math.sqrt(eig[0])),self.scaleRho(math.sqrt(eig[1])),0,360, 360*math.acos(Rinv[0][0])/math.pi)


    def drawJets(self, eV, color, lineWidth) :
        if not self.jets : return
        self.legendFunc(color, name = "cleanJet".join(self.jets), desc = "clean jets (%s%s)"%self.jets)
        
        jets = eV["CorrectedP4".join(self.jets)]
        fPU = eV["PileUpPtFraction".join(self.jets)]
        
        for iJet in eV["Indices".join(self.jets)] :
            jet = jets.at(iJet)
            icolor = r.kGray if fPU[iJet] > 0.7 else color
            self.rhoPhiPad.cd(); self.drawP4(self.rhoPhiCoords, jet, icolor, lineWidth, self.arrow.GetDefaultArrowSize() )
            self.etaPhiPad.cd(); self.drawCircle(jet, icolor, lineWidth, circleRadius = self.jetRadius)
            if jet.pt()>35 :
                self.drawCircle(jet, icolor, lineWidth, circleRadius = self.jetRadius*(1 + 0.0005*(jet.pt()-30)), lineStyle=1)
                self.drawCircle(jet, icolor, lineWidth, circleRadius = self.jetRadius*(1 - 0.0005*(jet.pt()-30)), lineStyle=1)
            
                    
    def drawLeptons(self, eV, color, lineWidth, kind = "muon", desc = "") :
        lepton = getattr(self,kind+'s')
        if not lepton : return
        self.legendFunc(color, name = kind.join(lepton) if not desc else desc, desc = kind + "s (%s%s)"%lepton if not desc else desc)
        p4 = eV["P4".join(lepton)]
        self.rhoPhiPad.cd()
        for i in eV["IndicesAnyIso".join(lepton)] :
            self.drawP4(self.rhoPhiCoords, p4.at(i), color, lineWidth, self.arrow.GetDefaultArrowSize() )
                        
    def preparePads(self) :
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


    def drawLegend(self) :
        self.legendPad.cd();
        
        legend = r.TLegend(0.0, 0.0, 1.0, 1.0)
        for item in self.legendList :
            self.line.SetLineColor(item[0])
            self.line.SetLineWidth(10)
            someArrow = self.line.DrawLine(0.0,0.0,0.0,0.0)
            legend.AddEntry(someArrow, item[1], item[2])
        legend.Draw("same")
        self.canvas.cd()
        self.legendPad.Draw()
        self.keep.append(legend)
        self.line.SetLineWidth(1)

    def legendFunc(self, color, name, desc) :
        if not self.legendDict[name] :
            self.legendDict[name] = True
            self.legendList.append( (color, self.renamedDesc(desc), "l") )

    def renamedDesc(self, desc) :
        if not self.prettyMode : return desc
        elif desc in self.prettyReName : return self.prettyReName[desc]
        else : return desc
        

    def prepareText(self, params, coords) :
        self.text.SetTextSize(params["size"])
        self.text.SetTextFont(params["font"])
        self.text.SetTextColor(params["color"])
        self.textSlope = params["slope"]

        self.textX = coords["x"]
        self.textY = coords["y"]
        self.textCounter = 0

    def printText(self, message, color = None) :
        if color != None : self.text.SetTextColor(color)
        self.text.DrawText(self.textX, self.textY - self.textCounter * self.textSlope, message)
        self.textCounter += 1

    def printEvent(self, eV, params, coords) :
        self.prepareText(params, coords)
        for message in ["Run   %#10d"%eV["run"],
                        "Ls    %#10d"%eV["lumiSection"],
                        "Event %#10d"%eV["event"],
                        "genQ (GeV) %#5.1f"%eV["genQ"] if not eV["isRealData"] else "",
                        ] :
            if message : self.printText(message)
        
    def printDiagnosis(self, eV, params, coords) :
        mode = eV["kinfitFailureModes"]
        self.prepareText(params, coords); self.printText('t    ' if 't' in mode else '', r.kGreen if 't' in mode and mode['t'] else r.kRed)
        self.prepareText(params, coords); self.printText('   /t' if '/t' in mode else '', r.kGreen if '/t' in mode and mode['/t'] else r.kRed)
        coords.update({"y":coords['y']-1.5*params['slope']});
        self.prepareText(params, coords)
        self.printText(eV["ttDecayMode"], r.kGreen if eV['ttDecayMode'] in ['ej','mj'] else r.kRed)
        for item in ['met','jet','','blep','nu','','had','bhad','glu'] :
            if item not in mode:
                coords.update({"x":coords['x']+0.27})
                self.prepareText(params, coords)
                continue
            self.printText(item if item in mode else '', r.kGreen if item in mode and mode[item] else r.kRed)
        
        


    def printNarrowText(self, eV) :
        self.narrowPad.Clear(); self.narrowPad.cd()
        defaults = {"size":0.085, "font":80, "color":r.kBlack, "slope":0.026}
        self.printEvent(   eV, params = defaults, coords = {"x":0.01, "y":0.98})
        if not eV['isRealData']: self.printDiagnosis(   eV, params = defaults, coords = {"x":0.01, "y":0.85})
        self.canvas.cd()
        self.narrowPad.Draw()




