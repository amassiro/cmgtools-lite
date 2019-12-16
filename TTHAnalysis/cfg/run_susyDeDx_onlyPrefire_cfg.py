##########################################################
##       CONFIGURATION FOR SUSY STOP SOFT B TREES       ##
##########################################################
#
# Prefire from https://github.com/CERN-PH-CMG/cmgtools-lite/blob/b3e71821ecf4ec50331f3b00d3109119d99d15c6/TTHAnalysis/python/analyzers/susyCore_modules_cff.py
#
# more properly: https://github.com/CERN-PH-CMG/cmgtools-lite/blob/94X_dev/TTHAnalysis/python/analyzers/PrefiringAnalyzer.py
#
#



import PhysicsTools.HeppyCore.framework.config as cfg
import re, sys

#-------- LOAD ALL ANALYZERS -----------
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption
from CMGTools.RootTools.samples.configTools import *
from CMGTools.RootTools.samples.autoAAAconfig import *


#
# Prefire
#

from CMGTools.TTHAnalysis.analyzers.PrefiringAnalyzer import PrefiringAnalyzer
PrefiringAnalyzer = cfg.Analyzer(
  PrefiringAnalyzer, name='PrefiringAnalyzer',
  #class_object= PrefiringAnalyzer,
  L1Maps = '$CMSSW_BASE/src/CMGTools/RootTools/data/L1PrefiringMaps_new.root',
  DataEra = '2017BtoF',
  UseJetEMPt = False ,
  PrefiringRateSystematicUncty =  0.2 , 
  SkipWarnings= True,
  photons='slimmedPhotons'
  )
  


from PhysicsTools.Heppy.analyzers.objects.all import *
## Photon Analyzer (generic)
photonAna = cfg.Analyzer(
    PhotonAnalyzer, name='photonAnalyzer',
    photons='slimmedPhotons',
    doPhotonScaleCorrections=False, 
    ptMin = 15,
    etaMax = 2.5,
    gammaID = "POG_PHYS14_25ns_Loose",
    rhoPhoton = 'fixedGridRhoFastjetAll',
    gamma_isoCorr = 'rhoArea',
    conversionSafe_eleVeto = False,
    do_mc_match = False,  # AM True
    do_randomCone = False,
)



# Select a list of good primary vertices (generic)
vertexAna = cfg.Analyzer(
    VertexAnalyzer, name="VertexAnalyzer",
    vertexWeight = None,
    fixedWeight = 1,
    verbose = False
    )


# Lepton Analyzer (generic)
#
# -> this is the producer of event.selectedLeptons 
#
lepAna = cfg.Analyzer(
    LeptonAnalyzer, name="leptonAnalyzer",
    # input collections
    muons='slimmedMuons',
    electrons='slimmedElectrons',
    rhoMuon= 'fixedGridRhoFastjetAll',
    rhoElectron = 'fixedGridRhoFastjetAll',
    # energy scale corrections and ghost muon suppression (off by default)
    doMuonScaleCorrections=False,
    doElectronScaleCorrections=False, # "embedded" in 5.18 for regression
    doSegmentBasedMuonCleaning=False,
    # inclusive very loose muon selection
    inclusive_muon_id  = "POG_ID_Loose",
    inclusive_muon_pt  = 3,
    inclusive_muon_eta = 2.4,
    inclusive_muon_dxy = 0.5,
    inclusive_muon_dz  = 1.0,
    muon_dxydz_track = "innerTrack",
    # loose muon selection
    loose_muon_id     = "POG_ID_Loose",
    loose_muon_pt     = 10,
    loose_muon_eta    = 2.4,
    loose_muon_dxy    = 0.05,
    loose_muon_dz     = 0.1,
    loose_muon_relIso = 1e9,
    loose_muon_isoCut     = lambda muon : muon.relIso03 < 0.25,
    # inclusive very loose electron selection
    inclusive_electron_id  = "",
    inclusive_electron_pt  = 5,
    inclusive_electron_eta = 2.5,
    inclusive_electron_dxy = 0.5,
    inclusive_electron_dz  = 1.0,
    inclusive_electron_lostHits = 1.0,
    # loose electron selection
    loose_electron_id     = "POG_Cuts_ID_FALL17_94X_v1_ConvVetoDxyDz_Veto",
    loose_electron_pt     = 10,
    loose_electron_eta    = 2.5,
    loose_electron_dxy    = 0.5,
    loose_electron_dz     = 1.0,
    loose_electron_relIso = 1e9,
    loose_electron_isoCut = lambda elec : elec.relIso03 < 0.25,
    loose_electron_lostHits = 10.0,
    # muon isolation correction method (can be "rhoArea" or "deltaBeta")
    mu_isoCorr = "deltaBeta" ,
    mu_effectiveAreas = "Fall17", #(can be 'Data2012' or 'Phys14_25ns_v1' or 'Spring15_25ns_v1')
    # electron isolation correction method (can be "rhoArea" or "deltaBeta")
    ele_isoCorr = "rhoArea" ,
    ele_effectiveAreas = "Fall17" , #(can be 'Data2012' or 'Phys14_25ns_v1' or 'Spring15_25ns_v1')
    ele_tightId = "Cuts_FALL17_94X_v1_ConvVetoDxyDz" ,
    # Mini-isolation, with pT dependent cone: will fill in the miniRelIso, miniRelIsoCharged, miniRelIsoNeutral variables of the leptons (see https://indico.cern.ch/event/368826/ )
    doMiniIsolation = False, # off by default since it requires access to all PFCandidates 
    packedCandidates = 'packedPFCandidates',
    miniIsolationPUCorr = 'rhoArea', # Allowed options: 'rhoArea' (EAs for 03 cone scaled by R^2), 'deltaBeta', 'raw' (uncorrected), 'weights' (delta beta weights; not validated)
    miniIsolationVetoLeptons = None, # use 'inclusive' to veto inclusive leptons and their footprint in all isolation cones
    doDirectionalIsolation = [], # calculate directional isolation with leptons (works only with doMiniIsolation, pass list of cone sizes)
    doFixedConeIsoWithMiniIsoVeto = False, # calculate fixed cone isolations with the same vetoes used for miniIso,
    # minimum deltaR between a loose electron and a loose muon (on overlaps, discard the electron)
    min_dr_electron_muon = 0.05,
    # do MC matching 
    do_mc_match = False, # AM True, # note: it will in any case try it only on MC, not on data
    do_mc_match_photons = "all",
    match_inclusiveLeptons = False, # match to all inclusive leptons
    )




## Jets Analyzer (generic)
jetAna = cfg.Analyzer(
    JetAnalyzer, name='jetAnalyzer',
    jetCol = 'slimmedJets',
    copyJetsByValue = False,      #Whether or not to copy the input jets or to work with references (should be 'True' if JetAnalyzer is run more than once)
    genJetCol = 'slimmedGenJets',
    rho = ('fixedGridRhoFastjetAll','',''),
    jetPt = 25.,
    jetEta = 4.7,
    jetEtaCentral = 2.4,
    cleanJetsFromLeptons = True,  
    jetLepDR = 0.4,
    jetLepArbitration = (lambda jet,lepton : lepton), # you can decide which to keep in case of overlaps; e.g. if the jet is b-tagged you might want to keep the jet
    cleanSelectedLeptons = True, #Whether to clean 'selectedLeptons' after disambiguation. Treat with care (= 'False') if running Jetanalyzer more than once
    minLepPt = 10,
    lepSelCut = lambda lep : True,
    relaxJetId = False,  
    doPuId = False, # Not commissioned in 7.0.X
    recalibrateJets = True, #'MC', # True, False, 'MC', 'Data'
    applyL2L3Residual = True, # Switch to 'Data' when they will become available for Data
    recalibrationType = "AK4PFchs",
    mcGT     = "Fall17_17Nov2017_V6_MC",
    dataGT   = [(1,"Fall17_17Nov2017B_V6_DATA"),(299337,"Fall17_17Nov2017C_V6_DATA"),(302030,"Fall17_17Nov2017D_V6_DATA"),(303435,"Fall17_17Nov2017E_V6_DATA"),(304911,"Fall17_17Nov2017F_V6_DATA")],
    jecPath = "${CMSSW_BASE}/src/CMGTools/RootTools/data/jec/",
    shiftJEC = 0, # set to +1 or -1 to apply +/-1 sigma shift to the nominal jet energies
    addJECShifts = False, # if true, add  "corr", "corrJECUp", and "corrJECDown" for each jet (requires uncertainties to be available!)
    jetPtOrUpOrDnSelection = False, # if true, apply pt cut on the maximum among central, JECUp and JECDown values of corrected pt
    smearJets = False,
    shiftJER = 0, # set to +1 or -1 to get +/-1 sigma shifts  
    alwaysCleanPhotons = False,
    cleanGenJetsFromPhoton = False,
    cleanJetsFromFirstPhoton = False,
    cleanJetsFromTaus = False,
    cleanJetsFromIsoTracks = False,
    doQG = False,
    do_mc_match = False,  # AM True
    collectionPostFix = "",
    calculateSeparateCorrections = True, # should be True if recalibrateJets is True, otherwise L1s will be inconsistent
    calculateType1METCorrection  = False,
    type1METParams = { 'jetPtThreshold':15., 'skipEMfractionThreshold':0.9, 'skipMuons':True },
    storeLowPtJets = False,
    )



from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer import *  

## Tree Producer
prefireTreeProducer = cfg.Analyzer(
    AutoFillTreeProducer, name='treeProducerPrefire',
    vectorTree = True, saveTLorentzVectors = False,  defaultFloatType = 'F', PDFWeights = [],
    globalVariables = [
        NTupleVariable("prefiringweight",      lambda ev: ev.prefiringweight, float, help="prefire weight"),
        NTupleVariable("prefiringweightup",    lambda ev: ev.prefiringweightup, float, help="prefire weight up"),
        NTupleVariable("prefiringweightdown",  lambda ev: ev.prefiringweightdown, float, help="prefire weight down"),
        ],
   )




# Core sequence of all common modules
coreSequence = [
  vertexAna,
  photonAna,
  lepAna,
  jetAna,
  PrefiringAnalyzer,
  prefireTreeProducer
]

#
#
#
  
## Configuration that depends on the region
region = getHeppyOption("region","sr") # use 'sr', 'cr1l'

## Sample production and setup
### Trigger
from CMGTools.RootTools.samples.triggers_13TeV_DATA2017 import *


### MC
from CMGTools.RootTools.samples.samples_13TeV_RunIIFall17MiniAODv2 import *

Wino_M_300_cTau_3  = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_3",  "Wino_M_300_cTau_3",  "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_300_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_10", "Wino_M_300_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_300_cTau_30 = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_30", "Wino_M_300_cTau_30", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_500_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_500_cTau_10", "Wino_M_500_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_650_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_650_cTau_10", "Wino_M_650_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_800_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_800_cTau_10", "Wino_M_800_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_1000_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_1000_cTau_10", "Wino_M_1000_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_500_cTau_20 = kreator.makeMCComponentFromEOS("Wino_M_500_cTau_20", "Wino_M_500_cTau_20", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_650_cTau_20 = kreator.makeMCComponentFromEOS("Wino_M_650_cTau_20", "Wino_M_650_cTau_20", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_800_cTau_20 = kreator.makeMCComponentFromEOS("Wino_M_800_cTau_20", "Wino_M_800_cTau_20", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_1000_cTau_20 = kreator.makeMCComponentFromEOS("Wino_M_1000_cTau_20", "Wino_M_1000_cTau_20", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)

Winos = [ Wino_M_300_cTau_3  , Wino_M_300_cTau_10 , Wino_M_300_cTau_30 , Wino_M_500_cTau_10, Wino_M_650_cTau_10, Wino_M_800_cTau_10, Wino_M_1000_cTau_10, Wino_M_500_cTau_20, Wino_M_650_cTau_20, Wino_M_800_cTau_20, Wino_M_1000_cTau_20 ]

Top = [ TTLep, TTHad, TTSemi] + Ts
VV  = [ WW, WZ, ZZ ]

Zll = [
    DYJetsM50_HT100to200,     DYJetsM50_HT100to200e,
    DYJetsM50_HT200to400,     DYJetsM50_HT200to400e,
    DYJetsM50_HT400to600,     DYJetsM50_HT400to600e,
    DYJetsM50_HT600to800,
    DYJetsM50_HT800to1200,
    #DYJetsM50_HT1200to2500,
    #DYJetsM50_HT2500toInf,
]

SelectedSamples = [
    #QCD_HT100to200,
    #TTSemi,
    #TBar_tch,
    WJets_HT100to200,
    WJets_HT600to800,
    ZvvJets_HT600to800,
]


if region == "sr":   
    mcSamples =  ( 
      #SelectedSamples
                   QCD
                 + Ws
                 + Zvv
                 + Zll 
                 + VV
                 + Top
                 )
    mcSignals = Winos
    mcTriggers = triggers_SOS_highMET[:] 
elif region == "cr1l": 
    mcSamples =  (
                   #QCD
                 #+ Ws
                 Zll 
                 #+ VV
                 #+ Top
                 )
    mcTriggers = triggers_1mu_iso + triggers_1e_iso + triggers_1e_noniso
    mcSignals = []

autoAAA(mcSamples)
cropToLumi(mcSamples, 10*41.7)
for c in mcSamples + mcSignals:
    c.triggers = mcTriggers


#
#
#


## Data
from CMGTools.RootTools.samples.samples_13TeV_DATA2017 import *
if region == "sr":   
    datasetsAndTriggers = [ ("MET", triggers_SOS_highMET) ]
elif region == "cr1l":   
    datasetsAndTriggers = [ ("SingleMuon", triggers_1mu_iso) #,
                          ]

json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
dataSamples = []; vetoTriggers = []
for (pdname, trigs) in datasetsAndTriggers:
    for d in dataSamples_31Mar2018:
        if pdname in d.name:
            d.json = json
            d.triggers = trigs[:]
            d.vetoTriggers = vetoTriggers[:]
            dataSamples.append(d)
    vetoTriggers += trigs

#
#
#


run = getHeppyOption("run","all")
if run == "all":    selectedComponents = mcSamples + mcSignals
elif run == "mc":   selectedComponents = mcSamples
elif run == "sig":  selectedComponents = mcSignals

if run == "mc":
    for c in  selectedComponents:
       c.splitFactor = len(c.files)/2


#-------- SEQUENCE -----------
sequence = cfg.Sequence( coreSequence )


#-------- HOW TO RUN -----------
test = getHeppyOption('test')
if test == "1":
    if getHeppyOption("sample"):
        selected = [ s for s in selectedComponents if getHeppyOption("sample") == s.name ]
        if not selected:
            print "Sample %s not found. Known samples are: %s" % [ getHeppyOption("sample"), ", ".join(sorted(s.name for s in selectedComponents)) ]
            sys.exit(1)
        selectedComponents = selected
    print "The test wil use %s " % selectedComponents[0].name
    selectedComponents = doTest1(selectedComponents[0], sequence=sequence, cache=True )
    print "The test wil use file %s " % selectedComponents[0].files[0]
elif test == "1E": 
    comp = dataSamples[0]
    #comp.files = [ ': /store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/90000/FB3B40F9-0C73-C144-BB08-8E2DFB7AD448.root ' ]
    comp.files = [ '/tmp/amassiro/FB3B40F9-0C73-C144-BB08-8E2DFB7AD448.root' ]
    selectedComponents = doTest1(comp, sequence=sequence, cache=False )
    comp.isMC = True
    print "The test wil use file %s " % comp.files[0]
    comp.triggers = [] 

printSummary(selectedComponents)
if getHeppyOption("justSummary",False): sys.exit(0)
config = autoConfig(selectedComponents, sequence)


