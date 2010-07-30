import ROOT as r
from analysisStep import analysisStep
#####################################
class jetPtSelector(analysisStep) :
    """jetPtSelector"""

    def __init__(self,jetCollection,jetSuffix,jetPtThreshold,jetIndex):
        self.jetPtThreshold = jetPtThreshold
        self.jetIndex = jetIndex
        self.jetP4s = jetCollection+"CorrectedP4"+jetSuffix
        
        self.moreName = "(%s; %s; corr. pT[%d]>=%.1f GeV)" % (jetCollection, jetSuffix, jetIndex, jetPtThreshold)

    def select (self,eventVars,extraVars) :
        p4s = eventVars[self.jetP4s]
        if p4s.size() <= self.jetIndex : return False
        return self.jetPtThreshold <= p4s.at(self.jetIndex).pt()
#####################################
class jetPtVetoer(analysisStep) :
    """jetPtVetoer"""

    def __init__(self,jetCollection,jetSuffix,jetPtThreshold,jetIndex):
        self.jetPtThreshold = jetPtThreshold
        self.jetIndex = jetIndex
        self.jetP4s = jetCollection+"CorrectedP4"+jetSuffix
        
        self.moreName = "(%s; %s; corr. pT[%d]<%.1f GeV)" % (jetCollection, jetSuffix, jetIndex, jetPtThreshold)

    def select (self,eventVars,extraVars) :
        p4s = eventVars[self.jetP4s]
        if p4s.size() <= self.jetIndex : return True
        return p4s.at(self.jetIndex).pt() < self.jetPtThreshold
#####################################
class leadingUnCorrJetPtSelector(analysisStep) :
    """leadingUnCorrJetPtSelector"""

    def __init__(self,jetCollectionsAndSuffixes,jetPtThreshold):
        self.jetCollectionsAndSuffixes = jetCollectionsAndSuffixes
        self.jetPtThreshold = jetPtThreshold

        self.moreName = "("+''.join(["%s%s;" % cS for cS in self.jetCollectionsAndSuffixes])
        self.moreName2 = "corr. pT[leading uncorr. jet]>=%.1f GeV" % self.jetPtThreshold 

    def select (self,eventVars,extraVars) :
        # Corrected pt of leading jet (by uncorrected pt) >= threshold
        for cS in self.jetCollectionsAndSuffixes :
            p4s = eventVars["%sCorrectedP4%s" % cS]
            corr = eventVars["%sCorrFactor%s" % cS]
            size = p4s.size()
            if not size : continue
            maxUncorrPt,index = max( [ (p4s.at(i).pt()/corr.at(i),i) for i in range(size) ] )
            if self.jetPtThreshold <= p4s.at(index).pt() :
                return True
        return False
#####################################
class cleanJetIndexProducer(analysisStep) :
    """cleanJetIndexProducer"""

    def __init__(self,jetCollection,jetSuffix,jetPtThreshold,jetEtaMax):
        self.jetCollection = jetCollection
        self.jetSuffix = jetSuffix
        self.jetPtThreshold = jetPtThreshold
        self.jetEtaMax = jetEtaMax
        
        self.moreName = "(%s; %s; jetID; " % (jetCollection, jetSuffix)
        self.moreName2 = " corr. pT>=%.1f GeV; |eta|<=%.1f)" % (jetPtThreshold,jetEtaMax)

        self.helper = r.cleanJetIndexHelper()
        self.helper.SetThresholds(self.jetPtThreshold,self.jetEtaMax)
        self.cleanJetIndices = r.std.vector('int')()
        self.otherJetIndices = r.std.vector('int')()
        self.cleanJetIndices.reserve(256)
        self.otherJetIndices.reserve(256)

    def uponAcceptance (self,eventVars,extraVars) :
        self.cleanJetIndices.clear()
        self.otherJetIndices.clear()

        p4Vector     =eventVars[self.jetCollection+'CorrectedP4'     +self.jetSuffix]
        emfVector    =eventVars[self.jetCollection+'EmEnergyFraction'+self.jetSuffix]
        fHpdVector   =eventVars[self.jetCollection+'JetIDFHPD'       +self.jetSuffix]
        n90HitsVector=eventVars[self.jetCollection+'JetIDN90Hits'    +self.jetSuffix]
        self.helper.Loop(p4Vector,emfVector,fHpdVector,n90HitsVector,self.cleanJetIndices,self.otherJetIndices)

        setattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix,[index for index in self.cleanJetIndices])
        setattr(extraVars,self.jetCollection+"otherJetIndices"+self.jetSuffix,[index for index in self.otherJetIndices])

        #setattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix,self.cleanJetIndices)
        #setattr(extraVars,self.jetCollection+"otherJetIndices"+self.jetSuffix,self.otherJetIndices)

#####################################
class cleanJetIndexProducerFromFlag(analysisStep) :
    """cleanJetIndexProducerFromFlag"""

    def __init__(self,jetCollection,jetSuffix,jetPtThreshold,jetEtaMax):
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.jetPtThreshold=jetPtThreshold
        self.jetEtaMax=jetEtaMax
        
        self.moreName="("+self.jetCollection+"; "+self.jetSuffix
        self.moreName+="; jetIDFlag; "

        self.moreName2=" "
        self.moreName2+="corr. "
        self.moreName2+="pT>="+str(self.jetPtThreshold)+" GeV"
        self.moreName2+="; |eta|<="+str(self.jetEtaMax)
        self.moreName2+=")"

        self.flagName="JetIDloose"
        if self.jetCollection[-2:]=="PF" : self.flagName="PF"+self.flagName
        
    def uponAcceptance (self,eventVars,extraVars) :
        p4Vector       =eventVars[self.jetCollection+'CorrectedP4'     +self.jetSuffix]
        jetIdFlagVector=eventVars[self.jetCollection+self.flagName     +self.jetSuffix]

        cleanString=self.jetCollection+"cleanJetIndices"+self.jetSuffix
        setattr(extraVars,cleanString,[])
        cleanJetIndices=getattr(extraVars,cleanString)

        otherString=self.jetCollection+"otherJetIndices"+self.jetSuffix
        setattr(extraVars,otherString,[])
        otherJetIndices=getattr(extraVars,otherString)

        size=p4Vector.size()
        
        for iJet in range(size) :
            #pt cut
            if p4Vector.at(iJet).pt()<self.jetPtThreshold : break #assumes sorted
            
            #if pass pt cut, add to "other" category
            otherJetIndices.append(iJet)
            
            #eta cut
            absEta=r.TMath.Abs(p4Vector.at(iJet).eta())
            if absEta>self.jetEtaMax : continue
            #jet ID cut
            if not jetIdFlagVector.at(iJet) : continue

            cleanJetIndices.append(iJet)
            otherJetIndices.remove(iJet)

        book = self.book(eventVars)
        print book
        book.fill(len(cleanJetIndices), cleanString, 15,-0.5,14.5, title=";number of jets passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin")
#####################################
class cleanJetEmfFilter(analysisStep) :
    """cleanJetEmfFilter"""

    def __init__(self,jetCollection,jetSuffix,jetPtThreshold,jetEmfMax):
        self.indicesName = "%scleanJetIndices%s" % jetCollection,jetSuffix
        self.p4sName = "%sCorrectedP4%s"         % jetCollection,jetSuffix
        self.emfName = "%sEmEnergyFraction%s"    % jetCollection,jetSuffix

        self.jetPtThreshold = jetPtThreshold
        self.jetEmfMax = jetEmfMax
        
        self.moreName = "(%s; %s" % self.jetCollection, self.jetSuffix
        self.moreName2 = " corr. pT>=%.1f GeV; EMF<=%.1f)" % jetPtThreshold,self.jetEmfMax

    def select (self,eventVars,extraVars) :
        cleanJetIndices = getattr( extraVars, self.indicesName)
        p4s = eventVars[self.p4sName]
        emf = eventVars[self.emfName]

        for index in cleanJetIndices :
            if p4s.at(index).pt() <= self.jetPtThreshold : #assumes sorted
                return True
            if emf.at(index) > self.jetEmfMax :
                return False
        return True
#####################################
class pfCleanJetIndexProducer(analysisStep) :
    """pfCleanJetIndexProducer"""

    def __init__(self,jetCollection,jetSuffix,jetPtThreshold,jetEtaMax):
        self.jetCollection = jetCollection
        self.jetSuffix = jetSuffix
        self.jetPtThreshold = jetPtThreshold
        self.jetEtaMax = jetEtaMax

        self.moreName = "(%s; %s; jetID; pT>=%.1f GeV; |eta|<=%.1f)" % jetCollection, jetSuffix, jetPtThreshold, jetEtaMax
        
    def step1 (self,eventVars,extraVars) :
        setattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix,[])
        self.cleanJetIndices = getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix)

        setattr(extraVars,self.jetCollection+"otherJetIndices"+self.jetSuffix,[])
        self.otherJetIndices = getattr(extraVars,self.jetCollection+"otherJetIndices"+self.jetSuffix)

        self.p4Vector = eventVars[self.jetCollection+'CorrectedP4'+self.jetSuffix]
        self.size = self.p4Vector.size()

    def jetLoop (self,eventVars) :
        for iJet in range(self.size) :
            #pt cut
            if self.p4Vector.at(iJet).pt() < self.jetPtThreshold : break #assumes sorted
            
            #if pass pt cut, add to "other" category
            self.otherJetIndices.append(iJet)
            
            #eta cut
            absEta = abs(self.p4Vector.at(iJet).eta())
            if self.jetEtaMax < absEta : continue
            
            self.cleanJetIndices.append(iJet)
            self.otherJetIndices.remove(iJet)

    def select (self,eventVars,extraVars) :
        self.step1(eventVars,extraVars)
        self.jetLoop(eventVars)
        return True
######################################
class nCleanJetHistogrammer(analysisStep) :
    """nCleanJetHistogrammer"""

    def __init__(self,jetCollection,jetSuffix):
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.moreName="("+self.jetCollection+" "+self.jetSuffix+")"

    def bookHistos(self) :
        nBins=15
        title=";number of jets passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin"
        self.nCleanJetsHisto=r.TH1D(self.jetCollection+"nCleanJets"+self.jetSuffix,title,nBins,-0.5,nBins-0.5)
        
    def uponAcceptance (self,eventVars,extraVars) :
        self.nCleanJetsHisto.Fill( len(getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix)) )
######################################
class nCleanJetEventFilter(analysisStep) :
    """nCleanJetEventFilter"""

    def __init__(self,jetCollection,jetSuffix,nCleanJets):
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.nCleanJets=nCleanJets
        self.moreName="("+self.jetCollection+" "+self.jetSuffix+">="+str(self.nCleanJets)+")"
        
    def select (self,eventVars,extraVars) :
        return len(getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix))>=self.nCleanJets
######################################
class nOtherJetEventFilter(analysisStep) :
    """nOtherJetEventFilter"""

    def __init__(self,jetCollection,jetSuffix,nOtherJets):
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.nOtherJets=nOtherJets
        self.moreName="("+self.jetCollection+" "+self.jetSuffix+"<"+str(self.nOtherJets)+")"
        
    def select (self,eventVars,extraVars) :
        return len(getattr(extraVars,self.jetCollection+"otherJetIndices"+self.jetSuffix))<self.nOtherJets
#####################################
class cleanJetHtMhtProducer(analysisStep) :
    """cleanJetHtMhtProducer"""

    def __init__(self,jetCollection,jetSuffix):
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.moreName="("
        self.moreName+=self.jetCollection
        self.moreName+="; "
        self.moreName+=self.jetSuffix
        self.moreName+=")"

        self.helper=r.htMhtHelper()
        self.cleanJetIndices=r.std.vector('int')()
        self.cleanJetIndices.reserve(256)
        
    def select (self,eventVars,extraVars) :
        p4Vector=eventVars[self.jetCollection+'CorrectedP4'+self.jetSuffix]
        cleanJetIndices=getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix)

        self.cleanJetIndices.clear()
        for index in cleanJetIndices :
            self.cleanJetIndices.push_back(index)

        self.helper.Loop(p4Vector,self.cleanJetIndices)
        setattr(extraVars,self.jetCollection+"Mht" +self.jetSuffix,self.helper.GetMht()  )
        setattr(extraVars,self.jetCollection+"Ht"  +self.jetSuffix,self.helper.GetHt()   )
        setattr(extraVars,self.jetCollection+"HtEt"+self.jetSuffix,self.helper.GetHtEt() )
        return True
#####################################
class cleanJetHtMhtHistogrammer(analysisStep) :
    """cleanJetHtMhtHistogrammer"""

    def __init__(self,jetCollection,jetSuffix,corrRatherThanUnCorr):
        self.histoMax=1.0e3
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.corrRatherThanUnCorr=corrRatherThanUnCorr
        self.corrString=""
        if (not self.corrRatherThanUnCorr) : self.corrString=" uncorrected"
        self.moreName="("
        self.moreName+=self.jetCollection
        self.moreName+="; "
        self.moreName+=self.jetSuffix
        self.moreName+=self.corrString
        self.moreName+=")"

    def bookHistos(self):
        self.ht_Histo          =r.TH1D(self.jetCollection+"ht"   +self.jetSuffix+" "+self.corrString
                                       ,";"+self.corrString+" H_{T} (GeV) from clean jet p_{T}'s;events / bin" ,50,0.0,self.histoMax)
        self.ht_et_Histo       =r.TH1D(self.jetCollection+"ht_et"+self.jetSuffix+" "+self.corrString
                                       ,";"+self.corrString+" H_{T} (GeV) from clean jet E_{T}'s;events / bin" ,50,0.0,self.histoMax)
        self.mht_Histo         =r.TH1D(self.jetCollection+"mht"  +self.jetSuffix+" "+self.corrString
                                       ,";"+self.corrString+" #slash{H}_{T} (GeV) from clean jets;events / bin",50,0.0,self.histoMax)
        self.m_Histo           =r.TH1D(self.jetCollection+"m"    +self.jetSuffix+" "+self.corrString
                                 ,";"+self.corrString+" mass (GeV) of system of clean jets;events / bin" ,50,0.0,self.histoMax)
        self.mHtOverHt_Histo   =r.TH1D(self.jetCollection+"mHtOverHt"   +self.jetSuffix+" "+self.corrString
                                 ,";"+self.corrString+" MHT / H_{T} (GeV) from clean jet p_{T}'s;events / bin" ,50,0.0,1.1)

        title=";"+self.corrString+" H_{T} (GeV) from clean jets;"+self.corrString+" #slash{H}_{T} (GeV) from clean jet p_{T}'s;events / bin"
        self.mhtHt_Histo=r.TH2D(self.jetCollection+"mht_vs_ht"+self.jetSuffix+" "+self.corrString
                                ,title,50,0.0,self.histoMax,50,0.0,self.histoMax)
        
    def uponAcceptance (self,eventVars,extraVars) :
        mhtVar="Mht"
        htVar="Ht"
        htEtVar="HtEt"
        
        if (not self.corrRatherThanUnCorr) :
            mhtVar+="UnCorr"
            htVar+="UnCorr"
            htEtVar+="UnCorr"
            
        mht =getattr(extraVars,self.jetCollection+mhtVar +self.jetSuffix)
        ht  =getattr(extraVars,self.jetCollection+htVar  +self.jetSuffix)
        htet=getattr(extraVars,self.jetCollection+htEtVar+self.jetSuffix)

        self.mht_Histo.Fill(mht.pt())
        self.ht_Histo.Fill(ht)
        self.ht_et_Histo.Fill(htet)
        self.m_Histo.Fill(mht.mass())
        self.mhtHt_Histo.Fill(ht,mht.pt())

        value=-1.0
        if (ht>0.0) : value=mht.pt()/ht
        self.mHtOverHt_Histo.Fill(value)
#####################################
class cleanJetPtHistogrammer(analysisStep) :
    """cleanJetPtHistogrammer"""

    def __init__(self,jetCollection,jetSuffix,fakeUnCorr=False) :
        self.histoMax=500.0
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.fakeUnCorr=fakeUnCorr
        self.moreName="("
        self.moreName+=self.jetCollection
        self.moreName+="; "
        self.moreName+=self.jetSuffix
        self.moreName+=")"

    def bookHistos(self) :
        self.ptAllHisto=          r.TH1D(self.jetCollection+" ptAll "          +self.jetSuffix
                                         ,";p_{T} (GeV) of clean jets;events / bin"                   ,50,0.0,self.histoMax)
        self.ptLeadingHisto=      r.TH1D(self.jetCollection+" ptLeading "      +self.jetSuffix
                                         ,";p_{T} (GeV) of leading clean jet;events / bin"            ,50,0.0,self.histoMax)
        if (not self.fakeUnCorr) :
            self.ptUnCorrAllHisto=    r.TH1D(self.jetCollection+" ptUnCorrAll "    +self.jetSuffix
                                             ,";uncorrected p_{T} (GeV) of clean jets;events / bin"       ,50,0.0,self.histoMax)
            self.ptUnCorrLeadingHisto=r.TH1D(self.jetCollection+" ptUnCorrLeading "+self.jetSuffix
                                             ,";uncorrected p_{T} (GeV) of leading clean jet;events / bin",50,0.0,self.histoMax)

    def uponAcceptance (self,eventVars,extraVars) :
        ptleading=0.0
        ptuncorrleading=0.0
        p4Vector=eventVars[self.jetCollection+'CorrectedP4'+self.jetSuffix]
        corrFactorVector=[1.0]*p4Vector.size()
        if not self.fakeUnCorr :
            corrFactorVector=eventVars[self.jetCollection+'CorrFactor'+self.jetSuffix]

        cleanJetIndices=getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix)

        for iJet in cleanJetIndices :
            p4=p4Vector.at(iJet)
            pt=p4.pt()
            pt_uc=pt
            if (not self.fakeUnCorr) : pt_uc/=corrFactorVector.at(iJet)

            self.ptAllHisto.Fill(pt)
            if (not self.fakeUnCorr) : self.ptUnCorrAllHisto.Fill(pt_uc)

            if (pt>ptleading) :
                ptleading=pt
            if (pt_uc>ptuncorrleading) :
                ptuncorrleading=pt_uc

        self.ptLeadingHisto.Fill(ptleading)
        if (not self.fakeUnCorr) : self.ptUnCorrLeadingHisto.Fill(ptuncorrleading)
#####################################
class cleanNJetAlphaProducer(analysisStep) :
    """cleanNJetAlphaProducer"""

    def __init__(self,jetCollection,jetSuffix):
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.moreName="("
        self.moreName+=self.jetCollection
        self.moreName+="; "
        self.moreName+=self.jetSuffix
        self.moreName+=")"

        self.helper=r.alphaHelper()
        self.cleanJetIndices=r.std.vector('int')()
        self.cleanJetIndices.reserve(256)

    def select (self,eventVars,extraVars) :
        nJetDeltaHt=0.0
        nJetAlphaT=0.0

        #return if fewer than two clean jets
        cleanJetIndices=getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix)
        if (len(cleanJetIndices)<2) :
            self.setExtraVars(extraVars,nJetDeltaHt,nJetAlphaT)
            return True

        #return if HT is tiny
        ht=getattr(extraVars,self.jetCollection+"Ht"+self.jetSuffix)
        if (ht<=1.0e-2) :
            self.setExtraVars(extraVars,nJetDeltaHt,nJetAlphaT)
            return True

        p4Vector=eventVars[self.jetCollection+"CorrectedP4"+self.jetSuffix]

        self.cleanJetIndices.clear()
        for index in cleanJetIndices :
            self.cleanJetIndices.push_back(index)

        self.helper.go(p4Vector,self.cleanJetIndices)
        nJetDeltaHt=self.helper.GetMinDiff()
        mht=getattr(extraVars,self.jetCollection+"Mht"+self.jetSuffix)        
        nJetAlphaT=0.5*(1.0-nJetDeltaHt/ht)/r.TMath.sqrt(1.0-(mht.pt()/ht)**2)

        self.setExtraVars(extraVars,nJetDeltaHt,nJetAlphaT)
        return True

    def setExtraVars(self,extraVars,nJetDeltaHt,nJetAlphaT) :
        setattr(extraVars,self.jetCollection+"nJetDeltaHt"+self.jetSuffix,nJetDeltaHt)
        setattr(extraVars,self.jetCollection+"nJetAlphaT"+self.jetSuffix,nJetAlphaT)
#####################################
class cleanDiJetAlphaProducer(analysisStep) :
    """cleanDiJetAlphaProducer"""

    def __init__(self,jetCollection,jetSuffix):
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.moreName="("
        self.moreName+=self.jetCollection
        self.moreName+="; "
        self.moreName+=self.jetSuffix
        self.moreName+=")"
        self.lvSum=r.Math.LorentzVector(r.Math.PxPyPzE4D('double'))(0.0,0.0,0.0,0.0)
        
    def select (self,eventVars,extraVars) :
        self.lvSum.SetCoordinates(0.0,0.0,0.0,0.0)

        diJetM       =0.0
        diJetMinPt   =1.0e6
        diJetMinEt   =1.0e6
        diJetAlpha   =0.0
        diJetAlpha_Et=0.0

        cleanJetIndices=getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix)
        #return if not dijet
        if (len(cleanJetIndices)!=2) :
            self.setExtraVars(extraVars,diJetM,diJetMinPt,diJetMinEt,diJetAlpha,diJetAlpha_Et)
            return True
            
        p4Vector=eventVars[self.jetCollection+"CorrectedP4"+self.jetSuffix]
        for iJet in cleanJetIndices :
            jet=p4Vector.at(iJet)
            pt=jet.pt()
            Et=jet.Et()

            if (pt<diJetMinPt) : diJetMinPt=pt
            if (Et<diJetMinEt) : diJetMinEt=Et

            self.lvSum+=jet

        diJetM=self.lvSum.M()
        
        if (diJetM>0.0) :
            diJetAlpha   =diJetMinPt/diJetM
            diJetAlpha_Et=diJetMinEt/diJetM

        self.setExtraVars(extraVars,diJetM,diJetMinPt,diJetMinEt,diJetAlpha,diJetAlpha_Et)
        return True

    def setExtraVars(self,extraVars,diJetM,diJetMinPt,diJetMinEt,diJetAlpha,diJetAlpha_Et) :
            setattr(extraVars,self.jetCollection+"diJetM"       +self.jetSuffix,diJetM)
            setattr(extraVars,self.jetCollection+"diJetMinPt"   +self.jetSuffix,diJetMinPt)
            setattr(extraVars,self.jetCollection+"diJetMinEt"   +self.jetSuffix,diJetMinEt)
            setattr(extraVars,self.jetCollection+"diJetAlpha"   +self.jetSuffix,diJetAlpha)
            setattr(extraVars,self.jetCollection+"diJetAlpha_Et"+self.jetSuffix,diJetAlpha_Et)
#####################################
class alphaHistogrammer(analysisStep) :
    """alphaHistogrammer"""

    def __init__(self,jetCollection,jetSuffix) :
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
        self.moreName="("
        self.moreName+=self.jetCollection
        self.moreName+="; "
        self.moreName+=self.jetSuffix
        self.moreName+=")"
        
    def bookHistos(self) :
        bins=100
        min=0.0
        max=2.0
        self.diJetAlpha_Histo   =r.TH1D(self.jetCollection+"dijet alpha"   +self.jetSuffix,";di-jet #alpha (using p_{T});events / bin"   ,bins,min,max)
        #self.diJetAlpha_ET_Histo=r.TH1D(self.jetCollection+"dijet alpha_ET"+self.jetSuffix,";di-jet #alpha (using E_{T});events / bin"   ,bins,min,max)
        self.nJetAlphaT_Histo   =r.TH1D(self.jetCollection+"njet alphaT"   +self.jetSuffix,";N-jet #alpha_{T} (using p_{T});events / bin",bins,min,max)
        self.nJetDeltaHt_Histo  =r.TH1D(self.jetCollection+"njet deltaHt"  +self.jetSuffix,";N-jet #Delta H_{T} (GeV);events / bin",50,0.0,500.0)
        self.alpha2D_c_Histo=r.TH2D("deltaHtOverHt vs mHtOverHt",
                                    ";#slash(H_{T}) / H_{T};#Delta H_{T} of two pseudo-jets / H_{T};events / bin",
                                    30,0.0,1.0,
                                    30,0.0,0.7)
        
    def uponAcceptance (self,eventVars,extraVars) :
        self.diJetAlpha_Histo.Fill(   getattr(extraVars,self.jetCollection+"diJetAlpha"   +self.jetSuffix))
        #self.diJetAlpha_ET_Histo.Fill(getattr(extraVars,self.jetCollection+"diJetAlpha_Et"+self.jetSuffix))
        self.nJetAlphaT_Histo.Fill(   getattr(extraVars,self.jetCollection+"nJetAlphaT"   +self.jetSuffix))

        mht=getattr(extraVars,self.jetCollection+"Mht"+self.jetSuffix).pt()
        ht=getattr(extraVars,self.jetCollection+"Ht"+self.jetSuffix)
        deltaHt=getattr(extraVars,self.jetCollection+"nJetDeltaHt"+self.jetSuffix)
        self.nJetDeltaHt_Histo.Fill(deltaHt)
        self.alpha2D_c_Histo.Fill(mht/ht,deltaHt/ht)
#####################################
class metHistogrammer(analysisStep) :
    """metHistogrammer"""

    def __init__(self,metCollection,tag) :
        self.tag=tag
        self.metCollection=metCollection
        
    def bookHistos(self) :
        bins=80
        min=0.0
        max=80.0
        self.caloMetNoHf_Histo=r.TH1D(self.metCollection+" "+self.tag,";"+self.metCollection+" p_{T} (GeV);events / bin",bins,min,max)
        
    def uponAcceptance (self,eventVars,extraVars) :
        self.caloMetNoHf_Histo.Fill(   eventVars[self.metCollection].pt() )
#####################################
class deltaPhiProducer(analysisStep) :
    """deltaPhiProducer"""

    def __init__(self,jetCollection,jetSuffix) :
        self.jetCollection=jetCollection
        self.jetSuffix=jetSuffix
    
    def select(self,eventVars,extraVars) :

        setattr(extraVars,self.jetCollection+"deltaPhi01"+self.jetSuffix, -4.0)
        setattr(extraVars,self.jetCollection+"deltaR01"  +self.jetSuffix,-40.0)
        setattr(extraVars,self.jetCollection+"deltaEta01"+self.jetSuffix,-40.0)

        p4Vector=eventVars[self.jetCollection+'CorrectedP4'+self.jetSuffix]
        cleanJetIndices=getattr(extraVars,self.jetCollection+"cleanJetIndices"+self.jetSuffix)

        if (len(cleanJetIndices)>=2) :
            jet0=p4Vector[cleanJetIndices[0]]
            jet1=p4Vector[cleanJetIndices[1]]
            setattr(extraVars,self.jetCollection+"deltaPhi01"+self.jetSuffix,r.Math.VectorUtil.DeltaPhi(jet0,jet1))
            setattr(extraVars,self.jetCollection+"deltaR01"  +self.jetSuffix,r.Math.VectorUtil.DeltaR(jet0,jet1))
            setattr(extraVars,self.jetCollection+"deltaEta01"+self.jetSuffix,jet0.eta()-jet1.eta())
        return True
#####################################
class deltaPhiSelector(analysisStep) :
    """deltaPhiSelector"""

    def __init__(self,jetCollection,jetSuffix,minAbs,maxAbs) :
        self.jetCollection = jetCollection
        self.jetSuffix = jetSuffix
        self.minAbs = minAbs
        self.maxAbs = maxAbs
    
    def select(self,eventVars,extraVars) :
        value = abs( getattr(extraVars,self.jetCollection+"deltaPhi01"+self.jetSuffix) )
        if (value<self.minAbs or value>self.maxAbs) : return False
        return True
#####################################
class mHtOverHtSelector(analysisStep) :
    """mHtOverHtSelector"""

    def __init__(self,jetCollection,jetSuffix,min,max) :
        self.jetCollection = jetCollection
        self.jetSuffix = jetSuffix
        self.min = min
        self.max = max
    
    def select(self,eventVars,extraVars) :
        mht = getattr(extraVars,self.jetCollection+"Mht"+self.jetSuffix).pt()
        ht = getattr(extraVars,self.jetCollection+"Ht"+self.jetSuffix)
        if (ht<1.0e-2) : return False
        value = mht/ht
        if (value<self.min or value>self.max) : return False
        return True
#####################################
class deltaPhiHistogrammer(analysisStep) :
    """deltaPhiHistogrammer"""

    def __init__(self,jetCollection,jetSuffix) :
        self.jetCollection = jetCollection
        self.jetSuffix = jetSuffix

    def bookHistos(self) :
        bins=50
        min=-4.0
        max= 4.0
        title=self.jetCollection+" deltaPhi01 "+self.jetSuffix
        self.deltaPhi01_Histo=r.TH1D(title,";"+title+";events / bin",bins,min,max)

        bins=20
        min= 0.0
        max=10.0
        title=self.jetCollection+" deltaR01 "+self.jetSuffix
        self.deltaR01_Histo=r.TH1D(title,";"+title+";events / bin",bins,min,max)

        bins=50
        min=-10.0
        max= 10.0
        title=self.jetCollection+" deltaEta01 "+self.jetSuffix
        self.deltaEta01_Histo=r.TH1D(title,";"+title+";events / bin",bins,min,max)
        
    def uponAcceptance (self,eventVars,extraVars) :
        self.deltaPhi01_Histo.Fill(getattr(extraVars,self.jetCollection+"deltaPhi01"+self.jetSuffix))
        self.deltaR01_Histo.Fill(getattr(extraVars,self.jetCollection+"deltaR01"+self.jetSuffix))
        self.deltaEta01_Histo.Fill(getattr(extraVars,self.jetCollection+"deltaEta01"+self.jetSuffix))
#####################################
