import collections, ROOT as r
from supy import utils,analysisStep
#####################################
try:
    import pdgLookup
    pdgLookupExists = True
except ImportError:
    pdgLookupExists = False
#####################################
class genJetPrinter(analysisStep) :

    def __init__(self, cs = None) :
        self.cs = cs
        self.moreName="(%s%s)"%self.cs
        self.p4 = "%sGenJetsP4%s"%self.cs

    def uponAcceptance (self,eventVars) :
        p4Vector = eventVars[self.p4]

        print " jet   pT (GeV)    eta    phi    emF   hadF   invF   auxF"
        print "---------------------------------------------------------"
        for iJet in range(len(p4Vector)) :
            jet=p4Vector[iJet]
            totalEnergy=jet.energy()
            
            outString  = " "
            outString += " %2d"   %iJet
            outString += "     %#6.1f"%jet.pt()
            outString += "   %#4.1f"%jet.eta()
            outString += "   %#4.1f"%jet.phi()
            print outString
        print
#####################################
class ParticleCountFilter(analysisStep) :
    def __init__(self, reqDict) :
        self.reqDict = reqDict
    def select (self,eventVars) :
        for key,value in self.reqDict.iteritems() :
            if eventVars["GenParticleCategoryCounts"][key]!=value : return False
        return True
#####################################
class scanHistogrammer(analysisStep) :
    def __init__(self, htVar = "", befOrAf = "", usePdfWeights = False) :
        self.tanBetaThreshold = 0.1
        for item in ["htVar"] :
            setattr(self, item, eval(item))
        self.moreName = self.htVar
        self.usePdfWeights = usePdfWeights
        self.befOrAf = befOrAf
#        self.m0Nbins = 81
#        self.m0Lo =      0.0
#        self.m0Hi =   2025.0
#
#        self.m12Nbins =  81
#        self.m12Lo =       0.0
#        self.m12Hi =    2025.0
#

## FOR T2CC, hardcoded :( ####
        self.m0Nbins = 34
        self.m0Lo =    90.0
        self.m0Hi =   260.0

        self.m12Nbins =  50
        self.m12Lo =     10.0
        self.m12Hi =    260.0


        self.bins = (self.m0Nbins, self.m12Nbins)
        self.lo = (self.m0Lo, self.m12Lo)
        self.hi = (self.m0Hi, self.m12Hi)

    def uponAcceptance (self, eventVars) :
        cat = eventVars["RA1category"]
        m0 = eventVars["susyScanmGL"]
        m12 = eventVars["susyScanmLSP"]
        xChi = "_%s"%eventVars["susyScanxCHI"] if ("susyScanxCHI" in eventVars) else ""
        title = ";m_{parent} (GeV);m_{LSP} (GeV)"
        pdfWeight =  1.0
        pdfSets = ["gencteq66","genMSTW2008nlo68cl","genNNPDF21","genct10"]
        if self.usePdfWeights :
            for pdfSet in pdfSets :
                if pdfSet in eventVars :
                	pdfWeights = eventVars[pdfSet]
                	for w,pdfWeight in enumerate(pdfWeights) :
                	    if self.befOrAf == "Before" : self.book.fill( (m0, m12), "nEvents_%s_%s"%(pdfSet,w), self.bins, self.lo, self.hi, pdfWeight, title = "%s;%s"%(title,""))
                	    else :
                	        self.book.fill( (m0, m12), "nEvents_"+cat+"_%s_%s"%(pdfSet,w), self.bins, self.lo, self.hi, pdfWeight, title = "%s;%s"%(title,""))
                	        self.book.fill( (m0, m12), "nEvents_%s_%s"%(pdfSet,w), self.bins, self.lo, self.hi, pdfWeight, title = "%s;%s"%(title,""))
        else :
            if self.befOrAf == "Before" : self.book.fill( (m0, m12), "nEvents%s"%(xChi),  self.bins, self.lo, self.hi, pdfWeight, title = "%s;%s"%(title,""))
            else :
                self.book.fill( (m0, m12), "nEvents_"+cat+"%s"%(xChi),  self.bins, self.lo, self.hi, pdfWeight, title = "%s;%s"%(title,""))
                self.book.fill( (m0, m12), "nEvents%s"%(xChi),  self.bins, self.lo, self.hi, pdfWeight, title = "%s;%s"%(title,""))
#####################################
class genParticleCountHistogrammer(analysisStep) :

    def __init__(self, tanBeta) :
        def nBins(lo, hi, stepSize) :
            return int(1+(hi-lo)/stepSize)

        self.tanBetaThreshold = 0.1
        self.tanBeta = tanBeta
        self.moreName = "tanBeta=%g"%self.tanBeta
        self.maxCountsPerCategory = 2 #0 ... this number counted explicitly; otherwise overflows

        #https://twiki.cern.ch/twiki/bin/view/CMS/SUSY38XSUSYScan#mSUGRA_Scans
        #Lo and Hi are both sampled in scan
        self.m0Lo =  10.0
        self.m0Hi = 500.0
        self.m0StepSize = 10.0
        self.nBinsM0 = nBins(self.m0Lo, self.m0Hi, self.m0StepSize)

        self.m12Lo = 100.0
        self.m12Hi = 350.0
        self.m12StepSize = 10.0
        self.nBinsM12 = nBins(self.m12Lo, self.m12Hi, self.m12StepSize)

        self.bins = (self.nBinsM0, self.nBinsM12)
        self.lo = (self.m0Lo-self.m0StepSize/2.0, self.m12Lo-self.m12StepSize/2.0)
        self.hi = (self.m0Hi+self.m0StepSize/2.0, self.m12Hi+self.m12StepSize/2.0)

        self.histoBaseName = "genParticleCounter"

    def makeCodeString(self,eventVars) :
        codeString = ""
        for category,count in eventVars["GenParticleCategoryCounts"].iteritems() :
            codeString += "_%s=%d"%(category, min(count, self.maxCountsPerCategory+1))
        return codeString
    
    def uponAcceptance (self,eventVars) :
        if abs(eventVars["susyScantanbeta"]-self.tanBeta)>self.tanBetaThreshold : return

        xs = eventVars["susyScanCrossSection"]
        m0 = eventVars["susyScanM0"]
        m12 = eventVars["susyScanM12"]
        #genMet = eventVars["metGenMetP4PF"].pt()
        codeString = self.makeCodeString(eventVars)        

        #self.book.fill( genMet, "GenMet", 40, 0.0, 800.0, title = ";gen. MET (GeV); events / bin")

        self.book.fill( (m0, m12), self.histoBaseName+codeString, self.bins, self.lo, self.hi,
                                   title = self.histoBaseName+codeString+";m_{0} (GeV);m_{1/2} (GeV)")

        #self.book.fill( (m0, m12), self.histoBaseName+codeString+"GenMet", self.bins, self.lo, self.hi,
        #                           w = genMet, title = self.histoBaseName+codeString+"GenMet;m_{0} (GeV);m_{1/2} (GeV)")
        #
        #self.book.fill( (m0, m12), self.histoBaseName+"GenMet", self.bins, self.lo, self.hi,
        #                           w = genMet, title = self.histoBaseName+"GenMet;m_{0} (GeV);m_{1/2} (GeV)")

        self.book.fill( (m0, m12), self.histoBaseName+"nEvents", self.bins, self.lo, self.hi,
                                   title = self.histoBaseName+"nEvents;m_{0} (GeV);m_{1/2} (GeV)")

        self.book.fill( (m0, m12), self.histoBaseName+"XS", self.bins, self.lo, self.hi,
                                   w = xs, title = self.histoBaseName+"XS;m_{0} (GeV);m_{1/2} (GeV)")
#####################################
class particlePrinter(analysisStep) :

    def __init__(self,minPt=-1.0,minStatus=-1):
        self.oneP4=utils.LorentzV()
        self.sumP4=utils.LorentzV()
        self.zeroP4=utils.LorentzV()
        self.minPt=minPt
        self.minStatus=minStatus
        
    def uponAcceptance (self,eventVars) :

        self.sumP4.SetCoordinates(0.0,0.0,0.0,0.0)

        mothers=set(eventVars["genMotherIndex"])
        print "pthat: ",eventVars["genpthat"]
        print "mothers: ",mothers
        print "-----------------------------------------------------------------------------------"
        print " i  st    mo         id            name        E        pt       eta    phi    mass"
        print "-----------------------------------------------------------------------------------"

        size=len(eventVars["genP4"])
        for iGen in range(size) :

            p4=eventVars["genP4"][iGen]
            if p4.pt()<self.minPt :
                continue

            status=eventVars["genStatus"][iGen]
            if status<self.minStatus :
                continue

            pdgId=eventVars["genPdgId"][iGen]
            outString=""
            outString+="%#2d"%iGen
            outString+=" %#3d"%status
            outString+="  %#4d"%eventVars["genMotherIndex"][iGen]
            outString+=" %#10d"%pdgId
            if pdgLookupExists : outString+=" "+pdgLookup.pdgid_to_name(pdgId).rjust(15)
            else :                 outString+="".rjust(16)
            outString+="  %#7.1f"%p4.E()
            outString+="  %#8.1f"%p4.pt()
            outString+="  %#8.1f"%p4.eta()
            outString+="  %#5.1f"%p4.phi()
            outString+="  %#6.1f"%p4.mass()
            #outString+="  %#5.1f"%p4.mass()
        
            if not (iGen in mothers) :
                outString+="   non-mo"
        #        self.sumP4+=self.oneP4
        #        #outString2="non-mo P4 sum".ljust(37)
        #        #outString2+="  %#7.1f"%self.sumP4.E()
        #        #outString2+="  %#8.1f"%self.sumP4.eta()
        #        #outString2+="  %#8.1f"%self.sumP4.pt()
        #        #outString2+="  %#5.1f"%self.sumP4.phi()
        #        #print outString2
        #
            print outString
        #
        #outString="non-mo P4 sum".ljust(37)
        #outString+="  %#7.1f"%self.sumP4.E()
        #outString+="  %#8.1f"%self.sumP4.eta()
        #outString+="  %#8.1f"%self.sumP4.pt()
        #outString+="  %#5.1f"%self.sumP4.phi()
        #print outString
        print
#####################################
#####################################
class topPrinter(analysisStep) :

    def uponAcceptance (self,ev) :
        ids = list(ev['genPdgId'])
        if not (6 in ids and -6 in ids) : return
        iTop = ids.index(6)
        iTbar= ids.index(-6)
        iQs = range(2,min(iTop,iTbar)) # interested in stuff listed between protons and tops
        mom = ev['genMotherIndex']
        p4s = ev['genP4']
        status = ev['genStatus']
        labels = {1:' d', -1:r'/d', 2:' u', -2:r'/u', 3:' s', -3:r'/s', 4:' c', -4:r'/c', 5:' b', -5:r'/b', 21:' g', 6:' t', -6:'/t'}

        iGs = filter(lambda i: ids[i]==21, range(max(iTop,iTbar)+1,len(ids)))
        iDaughters = filter(lambda i : mom[i]<min(iTop,iTbar), range(max(iTop,iTbar)+1,len(ids)))

        print '-'*50
        print '\t'.join(['item','mom','pt','m','eta','phi'])
        print
        for i in iQs+[iTop,iTbar]+iDaughters :
            print '\t'.join(str(d)[:5] for d in ["%d(%s)"%(i,labels[ids[i]] if ids[i] in labels else str(ids[i])), mom[i], p4s[i].pt(), '-', p4s[i].eta(), p4s[i].phi(), status[i]])
        print

        pD = sum([p4s[i] for i in iDaughters],utils.LorentzV())
        print '\t'.join(str(d)[:5] for d in ["Daught", '-', pD.pt(), pD.M(), pD.eta(), pD.phi()])
        p4 = p4s[2]+p4s[3]
        print '\t'.join(str(d)[:5] for d in ["[%s,%s]"%(2,3), '-', p4.pt(), p4.M(), p4.eta(), p4.phi()])
        p4 = p4s[4]+p4s[5]
        print '\t'.join(str(d)[:5] for d in ["[%s,%s]"%(4,5), '-', p4.pt(), p4.M(), p4.eta(), p4.phi()])
        tt = p4s[iTop] + p4s[iTbar]
        print '\t'.join(str(d)[:5] for d in ["ttbar", '-', tt.pt(), tt.M(), tt.eta(), tt.phi()])
        print
        print
        if abs(tt.E() - (p4s[4]+p4s[5]).E())>0.5 : print (50*' '), "2 -> 3+"
#####################################
class photonEfficiencyPlots(analysisStep) :

    def __init__(self, label, ptCut, etaCut, isoCut, deltaRCut, jets, photons) :
        for item in ["label","ptCut","etaCut","isoCut","deltaRCut","jets", "photons"] :
            setattr(self,item,eval(item))
        self.jetHt         = "%sSumEt%s"  %self.jets
        self.photonHt      = "%sSumEt%s"  %self.photons
        
        self.jetIndices    = "%sIndices%s"%self.jets
        self.photonIndices = "%sIndices%s"%self.photons

        self.moreName = ""
        if self.ptCut!=None     : self.moreName += "pT>%g GeV; "%self.ptCut
        if self.etaCut!=None    : self.moreName += "|eta|<%g; "%self.etaCut
        if self.isoCut!=None    : self.moreName += "iso<%g; "%self.isoCut
        if self.deltaRCut!=None : self.moreName += "deltaR>%g; "%self.deltaRCut

    def uponAcceptance (self, eventVars) :
        genP4s = eventVars["genP4"]
        nGen = genP4s.size()

        n = 0
        for genIndex in eventVars["genIndices"+self.label] :
            photon = genP4s.at(genIndex)
            pt = photon.pt()
            eta = photon.eta()
            phi = photon.phi()
            
            if pt<self.ptCut or self.etaCut<abs(eta) : continue

            if eventVars["category"+self.label][genIndex]=="otherMother" : continue

            deltaR = eventVars["genMinDeltaRPhotonOtherStatus3Photon"]
            if self.deltaRCut!=None and deltaR==None : continue
            if self.deltaRCut!=None and deltaR < self.deltaRCut : continue
            
            iso = eventVars["genIsolation"+self.label][genIndex]
            if self.isoCut!=None and self.isoCut < iso : continue

            n+=1
            self.book.fill(iso,"photonIso"+self.label, 100, 0.0,  100.0, title = ";gen photon isolation [5 GeV cut-off] (GeV);photons / bin")
            self.book.fill(eta,"photonEta"+self.label, 100, -3.0,   3.0, title = ";gen photon #eta;photons / bin")
            self.book.fill(pt, "photonPt"+self.label,  100,  0.0, 500.0, title = ";gen photon p_{T} (GeV);photons / bin")
            self.book.fill((eta, phi), "photonPhiVsEta"+self.label, (72, 72), (-3.0, -r.TMath.Pi()), (3.0, r.TMath.Pi()),
                                      title = ";gen photon #eta;gen photon #phi;photons / bin")

            nJets = len(eventVars[self.jetIndices])
            jetHt = eventVars[self.jetHt]

            photonIndices = eventVars[self.photonIndices]
            nPhotons      = len(photonIndices)
            photonHt      = eventVars[self.photonHt]
            
            self.book.fill(nJets,            "nJets"+self.label,              10, -0.5, 9.5,    title = ";nJets [gen photon satisfies cuts];photons / bin")
            self.book.fill(jetHt,            "jetHt"+self.label,             100,  0.0, 1000.0, title = ";H_{T} [jets] (GeV) [gen photon satisfies cuts];photons / bin")
            self.book.fill(nJets + nPhotons, "nJetsPlusnPhotons"+self.label,  10, -0.5, 9.5,    title = ";nJets+nPhotons [gen photon satisfies cuts];photons / bin")
            self.book.fill(jetHt + photonHt, "jetHtPlusPhotonHt"+self.label, 100,  0.0, 1000.0, title = ";H_{T} [jets+photons] (GeV) [gen photon satisfies cuts];photons / bin")
            if deltaR!=None : 
                self.book.fill(deltaR, "getMinDeltaRPhotonOtherStatus3Photon"+self.label,60, 0.0, 6.0,
                                          title = ";#DeltaR between st.3 photon and nearest daughterless st.3 particle; events / bin")
        self.book.fill(n,"nGenPhotons"+self.label, 10, -0.5, 9.5,title = ";N gen photons [gen photon satisfies cuts];photons / bin")
#####################################
class photonPurityPlots(analysisStep) :

    def __init__(self, label, jetCs, photonCs) :
        for item in ["label","jetCs","photonCs"] :
            setattr(self,item,eval(item))
        self.binLabels = ["photonMother", "quarkMother", "otherMother"]
        
    def uponAcceptance (self, eventVars) :
        genP4s   = eventVars["genP4"]
        jetSumP4 = eventVars["%sSumP4%s"%self.jetCs]

        photonIndices = eventVars["%sIndices%s"%self.photonCs]
        if not len(photonIndices) : return
        recoP4   = eventVars["%sP4%s"%self.photonCs].at( photonIndices[0] )
        categories = eventVars["category"+self.label]

        recoPt = recoP4.pt()
        matchedIndices = []
        for index in eventVars["genIndices"+self.label] :
            genP4 = genP4s.at(index)
            genPt  = genP4.pt()
            deltaR = r.Math.VectorUtil.DeltaR(recoP4,genP4)
            #if genPt>0.3*recoPt :
            #    self.book.fill(deltaR,"deltaRGenReco", 100, 0.0, 5.0, title = ";#DeltaR (GEN,RECO) photon when {gen pT > 0.3 reco pT};photons / bin")

            if deltaR>0.5 :
                continue
            matchedIndices.append(index)
            #self.book.fill(genPt,         categories[index]+"genPt" , 100, 0.0, 200.0, title = ";GEN photon p_{T} (GeV);photons / bin")
            #self.book.fill(recoPt,        categories[index]+"recoPt", 100, 0.0, 200.0, title = ";reco. photon p_{T} (GeV);photons / bin")
            #if jetSumP4 :
            #    self.book.fill(jetSumP4.pt(), categories[index]+"mht",    100, 0.0, 200.0, title = ";MHT (GeV);photons / bin")

        nMatch = len(matchedIndices)
        self.book.fill(nMatch,"nMatch",10, -0.5, 9.5, title = ";N gen. photons within #DeltaR 0.5 of the reco photon;photons / bin")
        label = "1"if len(matchedIndices)==1 else "gt1"
        for index in matchedIndices :
            self.book.fill( ( recoPt, genP4s.at(index).pt() ), "genVsRecoPt%s"%label,
                                       (50,50), (0.0,0.0), (200.0,200.0),
                                       title = "nMatch = %s;RECO photon p_{T} (GeV);GEN photon p_{T} (GeV);photons / bin"%label)
            self.book.fill(self.binLabels.index(categories[index]),"photonCategory%s"%label,
                                      len(self.binLabels), -0.5, len(self.binLabels)-0.5,
                                      title = ";photon category when nMatch = %s;photons / bin"%label)
#####################################
class genMotherHistogrammer(analysisStep) :

    def __init__(self, indexLabel, specialPtThreshold) :
        self.indexLabel = indexLabel
        self.specialPtThreshold = specialPtThreshold
        self.keyAll       = "motherIdVsPt%sAll"%self.indexLabel
        self.keyAllHighPt = "motherIdVsPt%sAllHighPt"%self.indexLabel        
        self.motherDict = collections.defaultdict(int)
        self.binLabels = []
        self.binLabels.append("other")

        self.addParticle( 1, "d"); self.addParticle(-1, "#bar{d}")
        self.addParticle( 2, "u"); self.addParticle(-2, "#bar{u}")
        self.addParticle( 3, "s"); self.addParticle(-3, "#bar{s}");
        self.addParticle( 4, "c"); self.addParticle(-4, "#bar{c}");
        self.addParticle(21, "gluon")
        self.addParticle(22, "photon")
        self.addParticle(111,"#pi^{0}")
        self.addParticle(221,"#eta")
        self.addParticle(223,"#omega")
        self.addParticle(331,"#eta^{/}")
        
    def addParticle(self, id, name) :
        self.binLabels.append(name)
        self.motherDict[id] = self.binLabels[-1]

    def fillSpecialHistos(self, eventVars, iParticle) :
        motherIndex = eventVars["genMother"].at(iParticle)
        #motherIndex = eventVars["genMotherIndex"].at(iParticle)
        p4 = eventVars["genP4"].at(iParticle)
        motherP4 = eventVars["genP4"].at(motherIndex)
        deltaRPhotonMother = r.Math.VectorUtil.DeltaR(p4,motherP4)
        deltaRPhotonOther  = r.Math.VectorUtil.DeltaR(p4,motherP4-p4)
        
        self.book.fill(motherP4.mass(), "mothersMass",
                                  20, -0.1, 0.4,
                                  title = ";mother's mass (GeV) [when GEN photon p_{T}> %.1f (GeV) and mother is u quark];photons / bin"%self.specialPtThreshold
                                  )
        self.book.fill(deltaRPhotonMother, "deltaRPhotonMother",
                                  20, 0.0, 1.5,
                                  title = ";#DeltaR(photon,mother) [when GEN photon p_{T}> %.1f (GeV) and mother is u quark];photons / bin"%self.specialPtThreshold
                                  )
        self.book.fill(deltaRPhotonOther, "deltaRPhotonOther",
                                  20, 0.0, 1.5,
                                  title = ";#DeltaR(photon,mother-photon) [when GEN photon p_{T}> %.1f (GeV) and mother is u quark];photons / bin"%self.specialPtThreshold
                                  )
        
    def uponAcceptance (self, eventVars) :
        indices = eventVars[self.indexLabel]
        if len(indices)==0 : return

        p4s = eventVars["genP4"]
        nBinsY = len(self.binLabels)
        for iParticle in indices :
            p4 = p4s.at(iParticle)
            pt = p4.pt()
            motherId = eventVars["genMotherPdgId"][iParticle]
            if not self.motherDict[motherId] :
                #print motherId,"not found"
                yValue = 0
            else :
                yValue = self.binLabels.index(self.motherDict[motherId])
            self.book.fill((pt,yValue), self.keyAll, (50,nBinsY), (0.0,-0.5), (500.0, nBinsY-0.5),
                                      title = ";GEN photon p_{T} (GeV);mother;photons / bin", yAxisLabels = self.binLabels)
            if pt>self.specialPtThreshold :
                self.book.fill(yValue, self.keyAllHighPt,
                                          nBinsY, -0.5, nBinsY-0.5,
                                          title = ";mother [when GEN photon p_{T}> %.1f (GeV)];photons / bin"%self.specialPtThreshold, xAxisLabels = self.binLabels)
                if motherId==2 : self.fillSpecialHistos(eventVars, iParticle)
#####################################
class zHistogrammer(analysisStep) :

    def __init__(self, jetCs) :
        self.jetCs = jetCs
        self.mhtName = "%sSumP4%s" % self.jetCs
        self.htName  = "%sSumEt%s"%self.jetCs
        
    def uponAcceptance (self,eventVars) :
        p4s = eventVars["genP4"]
        mht = eventVars[self.mhtName].pt()
        ht =  eventVars[self.htName]
        
        self.book.fill( len(cleanPhotonIndices), "photonMultiplicity", 10, -0.5, 9.5,
                        title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin"%self.cs)

        zIndices = eventVars["genIndicesZ"]
        if len(zS)>1 : return False
        for index in zIndices :
            Z = p4s.at(iPhoton)
            pt = photon.pt()


            self.book.fill(pt,           "%s%s%sPt" %(self.cs+(photonLabel,)), 50,  0.0, 500.0, title=";photon%s p_{T} (GeV);events / bin"%photonLabel)
            self.book.fill(photon.eta(), "%s%s%seta"%(self.cs+(photonLabel,)), 50, -5.0,   5.0, title=";photon%s #eta;events / bin"%photonLabel)

            self.book.fill((pt,mht), "%s%s%smhtVsPhotonPt"%(self.cs+(photonLabel,)),
                           (50, 50), (0.0, 0.0), (500.0, 500.0),
                           title=";photon%s p_{T} (GeV);MHT %s%s (GeV);events / bin"%(photonLabel,self.jetCs[0],self.jetCs[1])
                           )
            
            self.book.fill(mht/pt, "%s%s%smhtOverPhotonPt"%(self.cs+(photonLabel,)),
                           50, 0.0, 2.0, title=";MHT %s%s / photon%s p_{T};events / bin"%(self.jetCs[0],self.jetCs[1],photonLabel)
                           )

            #self.book.fill(pt-mht, "%s%s%sphotonPtMinusMht"%(self.cs+(photonLabel,)),
            #          100, -200.0, 200.0,
            #          title=";photon%s p_{T} - %s%sMHT (GeV);events / bin"%(photonLabel,self.jetCs[0],self.jetCs[1])
            #          )

            self.book.fill((pt-mht)/math.sqrt(ht+mht), "%s%s%sphotonPtMinusMhtOverMeff"%(self.cs+(photonLabel,)),
                           100, -20.0, 20.0,
                           title=";( photon%s p_{T} - MHT ) / sqrt( H_{T} + MHT )    [ %s%s ] ;events / bin"%(photonLabel,self.jetCs[0],self.jetCs[1])
                           )
#####################################
class DPhiHistogrammer(analysisStep) :
    
    def uponAcceptance(self,e) :
        mhtoverht = e["nonSusyFromSusySumP4PtOvernonSusyFromSusySumEt"]
        mht = e["nonSusyFromSusySumP4"].pt()
        ht = e["nonSusyFromSusySumEt"]
        dphiLsp = e["DPhiNeutralino"]
        p4 = e["genP4"]

        for iDaughter in  e["genIndicesNeutralino"] :
            iMother = e["genMotherIndex"].at(iDaughter)
            if iMother < 0 : continue
            dphi = abs(r.Math.VectorUtil.DeltaPhi(p4.at(iDaughter),p4.at(iMother)))
            self.book.fill(dphi, "DphiNeutralinoMother", 20, 0.0, r.TMath.Pi(),title=";#Delta#Phi between lsp and mother;lsp / bin")

            self.book.fill((dphi,mhtoverht), "MHT_HT_v_Dphilsp-mother",
                           (20, 50), (0.0, 0.0), (r.TMath.Pi(), 1.0),
                           title=";#Delta#Phi between LSP and mother;MHT/H_{T} from Gen Jets; LSP / bin"
                           )

#            self.book.fill((mhtoverht,dphi), "Dphilsp-mother_v_MHT_HT",
#                           (50,20), (0.0, 0.0), (1.0, r.TMath.Pi()),
#                           title=";MHT/H_{T} from Gen Jets;#Delta#Phi between LSP and mother; LSP / bin"
#                           )

            self.book.fill((dphi,mht), "MHT_v_DphiLSP-mother",
                           (20, 50), (0.0, 0.0), (r.TMath.Pi(), 1000),
                           title=";#Delta#Phi between LSP and mother;MHT (GeV) from Gen Jets; LSP / bin"
                           )

            self.book.fill((dphi,ht), "HT_v_DphiLSP-mother",
                           (20, 100), (0.0, 0.0), (r.TMath.Pi(), 2000),
                           title=";#Delta#Phi between LSP and mother;HT (GeV) from Gen Jets; LSP / bin"
                           )

            self.book.fill((dphiLsp,dphi), "DphiLSP-mother_v_DphiLsp",
                           (20, 20), (0.0, 0.0), (r.TMath.Pi(), r.TMath.Pi()),
                           title=";#Delta#Phi between LSP and mother;#Delta#Phi between LSPs; LSP / bin"
                           )

#####################################
class bqDotHistogrammer(analysisStep) :
    def f(self, x) :
        return x/5901.2 - 1.0

    def fill(self, x, y, sign = "") :
        self.book.fill((x, y), "bqDotHistogram%s"%sign, (100, 100), (-1.5, -1.5), (1.5, 1.5),
                       title = "W%s#rightarrowq#bar{q};(b.q)/5901.2 - 1.0;(b.#bar{q})/5901.2 - 1.0;events / bin"%sign)

    def uponAcceptance (self,eventVars) :
        p4s = eventVars["genP4"]
        for sign in ["+","-"] :
            b = p4s.at(eventVars["genIndicesb%s_t%sDaughters"%(sign,sign)][0])
            i0,i1 = eventVars["genIndicesW%sDaughters"%sign]
            if eventVars["genPdgId"][i0]>0 :
                q    = p4s.at(i0)
                qbar = p4s.at(i1)
            else :
                qbar = p4s.at(i0)
                q    = p4s.at(i1)
            self.fill(self.f(b.Dot(q)), self.f(b.Dot(qbar)), sign)
