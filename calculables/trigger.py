from supy import calculables,utils,wrappedChain
import collections,os, ROOT as r
import itertools,operator,random,math
try:
    import numpy as np
except:
    pass

#####################################
class lowestUnPrescaledTrigger(wrappedChain.calculable) :
    def __init__(self, sortedListOfPaths = []) :
        self.sortedListOfPaths = sortedListOfPaths
        self.cached = dict()
        self.moreName = "lowest unprescaled of " + utils.contract(self.sortedListOfPaths)

    def update(self, ignored) :
        key = (self.source["run"],self.source["lumiSection"])
        if key not in self.cached :
            self.cached[key] = None
            for path in self.sortedListOfPaths :
                if self.source["prescaled"][path]==1 :
                    self.cached[key] = path
                    break
        self.value = self.cached[key]
##############################
class TriggerWeight(calculables.secondary) :
    '''Probability that the event was triggered.'''
    
    def __init__(self, samples, tag = None, collection = ('muon','PF'), var = 'Pt', index = 'TriggeringIndex', triggers = [], thresholds = [], unreliable = {}) :
        self.fixes = collection
        self.var = var.join(collection)
        self.index = index.join(collection)
        self.binning = (200,0,200)
        for item in ['triggers','thresholds','tag','samples','unreliable'] : setattr(self,item,eval(item))

    def onlySamples(self) : return self.samples

    def setupDicts(self) :
        dd = collections.defaultdict
        self.epochLumis = dd(lambda : dd(set)) #set of lumis by run by (prescale set)
        self.lumis = dd(set) #set of lumis by run

    def setup(self,*_) :
        self.setupDicts()
        hists = self.fromCache( self.samples, ["susyTreeLumi",self.name], self.tag)
        if not all(all(hs.values()) for hs in hists.values()) :
            self.weights = r.TH1D('empty','',1,0,1)
            for i in range(3) : self.weights.SetBinContent(i, 1)
            return
        lumis = [ hists[sample]['susyTreeLumi'].GetBinContent(1) for sample in self.samples ]
        weights = [ hists[sample][self.name] for sample in self.samples ]
        for lumi,weight in zip(lumis,weights) : weight.Scale(lumi)
        self.weights = weights[0].Clone("random")
        for w in weights[1:] : self.weights.Add(w)
        self.weights.Scale(1./sum(lumis))
        for i in range(self.weights.GetNbinsX()+2) : assert 0 <= self.weights.GetBinContent(i) <= 1

    def triggerFired(self, val, triggered) :
        for thresh,trig in zip(self.thresholds,self.triggers) :
            if thresh > val : return False
            if triggered[trig] and (trig not in self.unreliable or
                                    self.source['prescaled'][trig] not in self.unreliable[trig]):
                return True
        return False

    def update(self,_) :
        index = self.source[self.index]
        self.value = (0 if index is None else
                      self.weights.GetBinContent(self.weights.FindFixBin(self.source[self.var][index])) if not self.source['isRealData'] else
                      1 if self.triggerFired(self.source[self.var][index], self.source['triggered']) else None )

    def uponAcceptance(self,ev) :
        index = ev[self.index]
        val = ev[self.var][index]
        
        self.book.fill(ev[self.name], 'value_%s'%self.name, 100, 0, 1, title = ';%s;events / bin'%self.name)
        self.book.fill(val, self.var,              *self.binning, title = ';%s;events / bin'%self.var)
        self.book.fill(val, self.var+'unweighted', *self.binning, title = ';%s;events / bin'%self.var, w = 1)

        if not ev['isRealData'] : return
        
        run,lumi = (ev["run"],ev["lumiSection"])
        if lumi in self.lumis[run] : return
        key = tuple([ev["prescaled"][trig] for trig in self.triggers])
        if self.unreliable : key = tuple([0 if trig in self.unreliable and ps in self.unreliable[trig] else ps for ps,trig in zip(key,self.triggers)])

        self.epochLumis[key][run].add(lumi)
        self.lumis[run].add(lumi)

    @property
    def staticEpochLumis(self) : return dict((epoch,dict(self.epochLumis[epoch].iteritems())) for epoch in self.epochLumis)

    def varsToPickle(self) : return ["staticEpochLumis"]
    def outputSuffix(self) : return "_%s/"%self.name

    def mergeFunc(self,products) :
        self.setupDicts()
        for eLumis in products["staticEpochLumis"] :
            for epoch in eLumis :
                for run in eLumis[epoch] :
                    self.epochLumis[epoch][run] |= eLumis[epoch][run]
        if not self.epochLumis : return

        lumiDir = self.outputFileName
        if not os.path.exists(lumiDir) : utils.mkdir(lumiDir)
        
        lumis = utils.luminosity.recordedInvMicrobarnsShotgun( [utils.jsonFromRunDict(self.epochLumis[epoch]) for epoch in self.epochLumis ] , cores = 4, cacheDir = lumiDir)
        probs = [ [1./prescale if prescale else 0 for prescale in epoch] for epoch in self.epochLumis ]
        inclu = [ [utils.unionProbability(prob[:i+1]) for i in range(len(prob))] for prob in probs ]
        thresholds = sorted(set(self.thresholds))
        inclusives = [ [max([0]+[p for p,t in zip(inc,self.thresholds) if t is thresh]) for thresh in thresholds] for inc in inclu ]
        weights = np.array(lumis).dot(inclusives) / sum(lumis)

        # write total lumi and weight by threshold
        lumiHist = r.TH1D('susyTreeLumi','luminosity from susyTree;;1/pb',1,0,1)
        lumiHist.SetBinContent(1,sum(lumis))
        lumiHist.Write()

        weightHist = r.TH1D(self.name,";%s;%s"%(self.var,self.name), len(thresholds), np.array([0.]+thresholds+[max(thresholds)+10]))
        for i,w in enumerate(weights) : weightHist.SetBinContent(i+2, w)
        weightHist.SetBinContent(len(thresholds)+3, 1)
        weightHist.Write()
############################################
class CrossTriggerWeight(calculables.secondary) :
    '''Implemented from description in CMS AN-2011/439 (4.6)'''

    def __init__(self, samples, triggers, jets = None) :
        self.jets = jets
        self.triggers = triggers
        self.samples = samples
        assert all([ abs(self.prob(effs)-probEffs) < 1e-15
                    for effs,probEffs in [(3*[1],1),
                                          (4*[1],1),
                                          (5*[1],1),
                                          (3*[0.5],0.5**3),
                                          ([0.95,0.95,0.1,0.1], sum( [0.95*0.95*0.1*0.1,
                                                                      0.95*0.95*0.1*(1-0.1),
                                                                      0.95*0.95*(1-0.1)*0.1,
                                                                      0.95*(1-0.95)*0.1*0.1,
                                                                      (1-0.95)*0.95*0.1*0.1] )),
                                          ([0.9,0.8,0.7,0.6,0.5], sum([0.9*0.8*0.7*0.6*0.5,
                                                                       0.9*0.8*0.7*0.6*(1-0.5),
                                                                       0.9*0.8*0.7*(1-0.6)*0.5,
                                                                       0.9*0.8*(1-0.7)*0.6*0.5,
                                                                       0.9*(1-0.8)*0.7*0.6*0.5,
                                                                       (1-0.9)*0.8*0.7*0.6*0.5,
                                                                       0.9*0.8*0.7*(1-0.6)*(1-0.5),
                                                                       0.9*0.8*(1-0.7)*0.6*(1-0.5),
                                                                       0.9*(1-0.8)*0.7*0.6*(1-0.5),
                                                                       (1-0.9)*0.8*0.7*0.6*(1-0.5),
                                                                       0.9*0.8*(1-0.7)*(1-0.6)*0.5,
                                                                       0.9*(1-0.8)*0.7*(1-0.6)*0.5,
                                                                       (1-0.9)*0.8*0.7*(1-0.6)*0.5,
                                                                       0.9*(1-0.8)*(1-0.7)*0.6*0.5,
                                                                       (1-0.9)*0.8*(1-0.7)*0.6*0.5,
                                                                       (1-0.9)*(1-0.8)*0.7*0.6*0.5
                                                                       ]))
                                          ]
                    ])
        
    def onlySamples(self) : return self.samples
    
    @staticmethod
    def gompertz((a,b,c),pT) : return a * math.exp( b * math.exp( c * pT ) )

    @staticmethod
    def abcs(triggertype) :
        return {"CentralJet30_A": ( (0.987,    -52.43,-0.166) ,
                                    (0.975,   -115.59,-0.187) ),
                "CentralJet30_B": ( (0.985,    -40.58,-0.160) ,
                                    (0.975,    -61.05,-0.163) ),
                "CentralPFJet30": ( (0.983,-189093.65,-0.410) ,
                                    (0.978,  -5914.00,-0.272) ),
                }[triggertype]

    def random_trigger_type(self) :
        rando = random.random()
        return ["CentralJet30_A",
                "CentralJet30_B",
                "CentralPFJet30"][ next( i for i,cumprob in enumerate(self.trigTypeCumProbs) if rando < cumprob) ]

    @staticmethod
    def prob(eff) :
        indices = range(len(eff))
        return sum( ( reduce(operator.mul, [eff[i] for i in iPass] + [1-eff[i] for i in indices if i not in iPass])
                    for njetsPass in range(3,len(eff)+1)
                    for iPass in itertools.combinations(indices, njetsPass) ) )

    def mcTriggeringProb(self) :
        indices = self.source['Indices'.join(self.jets)]
        if len(indices) < 3 : return 0
        p4 = self.source['CorrectedP4'.join(self.jets)]
        abcs = self.abcs( self.random_trigger_type() )
        return self.prob( [ self.gompertz( abcs[ 1.4 < abs(p4[i].eta()) ], p4[i].pt())
                            for i in indices] )

    def triggerFired(self) :
        return ( len(self.source['Indices'.join(self.jets)]) > 2 and 
                 any(self.source['triggered'][path] for path in self.triggers) )

    def update(self,_) :
        self.value = ( self.mcTriggeringProb() if not self.source['isRealData'] else
                       1                       if self.triggerFired() else
                       0 )

    def uponAcceptance(self,ev) :
        jets = self.source['CorrectedP4'.join(self.jets)]
        for i,iJet in list(enumerate(self.source['Indices'.join(self.jets)]))[:6] :
            barend = "barrel" if abs(jets[iJet].eta())<1.4 else 'endcap'
            pt = jets[iJet].pt()
            name = "jet%d_%s"%(i,barend)
            title = ";jet_{%d} pt (%s); events / bin"%(i,barend)
            self.book.fill( pt, name, 100,0,100, title = title )
            self.book.fill( pt, name+'_unweighted', 100,0,100,
                            title = 'unweighted '+title, w = 1)

        if not ev['isRealData'] : return
        trigger = next(ev['triggered'][path] for path in self.triggers)

        run,lumi = (ev["run"],ev["lumiSection"])
        if lumi in self.lumis[run] : return
        self.lumis[run].add(lumi)
        ps = ev['prescaled']
        self.prescales[(run,lumi)] = ( (1,0,0) if run < 174000 else
                                       ( 0,
                                         next(iter(sorted(ps[trig] for trig in self.triggers if 'PF' not in trig and ps[trig])), 0 ),
                                         next(iter(sorted(ps[trig] for trig in self.triggers if 'PF'     in trig and ps[trig])), 0 )) )

    def varsToPickle(self) : return ["prescales"]
    def outputSuffix(self) : return "_%s/"%self.name

    def mergeFunc(self,products) :
        dd = collections.defaultdict
        epochs = dd(lambda: dd(set))
        prescales = {}
        for ps in products['prescales'] : prescales.update(ps)
        if not prescales : return
        for (run,lumi),val in prescales.iteritems() : epochs[val][run].add(lumi)

        lumiDir = self.outputFileName
        if not os.path.exists(lumiDir) : utils.mkdir(lumiDir)

        lumis = utils.luminosity.recordedInvMicrobarnsShotgun( [ utils.jsonFromRunDict(rundict) for rundict in epochs.values() ],
                                                               cores = 4, cacheDir = lumiDir)
        lumiA,lumiB,lumiBPF = [ sum(lumis) for lumis in zip(*[(lumi/A if A else 0,
                                                               lumi/B if B else 0,
                                                               lumi/BPF if BPF else 0) for (A,B,BPF),lumi in zip(epochs.keys(),lumis) ] ) ]
        lumiHist = r.TH1D('lumis','luminosity by path;;1/pb', 3,0,3)
        lumiHist.SetBinContent(1,lumiA * 1e-6)   ; lumiHist.GetXaxis().SetBinLabel(1,'2011A')
        lumiHist.SetBinContent(2,lumiB * 1e-6)   ; lumiHist.GetXaxis().SetBinLabel(2,'2011B')
        lumiHist.SetBinContent(3,lumiBPF*1e-6) ; lumiHist.GetXaxis().SetBinLabel(3,'2011BPF')
        lumiHist.Write()
        return

    def setup(self,*_) :
        self.lumis = collections.defaultdict(set)
        self.prescales = {}

        hists = self.fromCache( self.samples, ["lumis"], tag = None )
        lumiHists = [shists['lumis'] for shists in hists.values()]
        if not all(lumiHists) :
            self.trigTypeCumProbs = 3*[1]
            return
        lumis = sum( [np.array([h.GetBinContent(i) for i in [1,2,3]]) for h in lumiHists ])
        self.trigTypeCumProbs = np.cumsum( lumis / sum(lumis) )
        return
