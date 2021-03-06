#!/bin/env python
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer import *  
from PhysicsTools.Heppy.analyzers.core.autovars import *  
from PhysicsTools.Heppy.analyzers.objects.autophobj import  *
from math import *

##------------------------------------------  
## LEPTON
##------------------------------------------  

leptonTypeXtracks = NTupleObjectType("leptonXtracks", baseObjectTypes = [ leptonType ], variables = [
    NTupleVariable("etaSc", lambda x : x.superCluster().eta() if abs(x.pdgId())==11 else -100, help="Electron supercluster pseudorapidity"),
    NTupleVariable("mcMatchPdgId",  lambda x : x.mcLep.pdgId() if getattr(x,'mcLep',None)!=None else -99, int, mcOnly=True, help="Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z): pdgId of the matched gen-level lepton, zero if non-prompt or fake"),
    NTupleVariable("mcPromptGamma", lambda x : x.mcPho.isPromptFinalState() if getattr(x,"mcPho",None) else 0, int, mcOnly=True, help="Photon isPromptFinalState"),
    NTupleVariable("triggerMatched", lambda x : 1 if x.matchedTrgObjSingleLep else 0, int, mcOnly=False, help="Is the lepton matched to a single lepton trigger"),
])
leptonTypeXtracks.removeVariable("relIsoAn04")
leptonTypeXtracks.removeVariable("eleCutIdSpring15_25ns_v1")
leptonTypeXtracks.removeVariable("ICHEPsoftMuonId")
leptonTypeXtracks.removeVariable("ICHEPmediumMuonId")

##------------------------------------------  
## TAU
##------------------------------------------  

tauTypeXtracks = NTupleObjectType("tauXtracks",  baseObjectTypes = [ tauType ], variables = [
     NTupleVariable("idMVAdR03", lambda x : x.idMVAdR03, int, help="1,2,3,4,5,6 if the tau passes the very loose to very very tight WP of the IsolationMVArun2v1DBdR03oldDMwLT discriminator"),
])
tauTypeXtracks.removeVariable("idMVANewDM")
tauTypeXtracks.removeVariable("idCI3hit")


jetTypeXtracks = NTupleObjectType("jetXtracks",  baseObjectTypes = [ jetTypeExtra ], variables = [
    NTupleVariable("chHEF", lambda x : x.chargedHadronEnergyFraction(), float, mcOnly = False, help="chargedHadronEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("neHEF", lambda x : x.neutralHadronEnergyFraction(), float, mcOnly = False,help="neutralHadronEnergyFraction (relative to uncorrected jet energy)"),
])
jetTypeXtracks.removeVariable("btagCMVA")
jetTypeXtracks.removeVariable("nLeptons")

jetTypeXtracksFwd = NTupleObjectType("jetFwd",  baseObjectTypes = [ jetType ], variables = [
    NTupleVariable("chHEF", lambda x : x.chargedHadronEnergyFraction(), float, mcOnly = False, help="chargedHadronEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("neHEF", lambda x : x.neutralHadronEnergyFraction(), float, mcOnly = False,help="neutralHadronEnergyFraction (relative to uncorrected jet energy)"),
])
jetTypeXtracks.removeVariable("btagCMVA")
jetTypeXtracksFwd.removeVariable("nLeptons")

##------------------------------------------  
## MET
##------------------------------------------  
  
metTypeXtracks = NTupleObjectType("metXtracks", baseObjectTypes = [ metType ], variables = [
])

metTypeXtracksBasic = NTupleObjectType("metXtracksBasic", baseObjectTypes = [ fourVectorType ], variables = [
])

##------------------------------------------  
## IsoTrackDeDx
##------------------------------------------

def _isoDBeta(x):
    return x.chargedHadronIso() + max(x.neutralHadronIso()+x.photonIso()-0.5*x.puChargedHadronIso(), 0)

isoTrackTypeDeDx = NTupleObjectType("isoTrackTypeDeDx", baseObjectTypes = [ particleType ], variables = [
    NTupleVariable("charge",   lambda x : x.charge(), int),
    NTupleVariable("dxy",   lambda x : x.dxy(), help="d_{xy} with respect to PV, in cm (with sign)"),
    NTupleVariable("dz",    lambda x : x.dz() , help="d_{z} with respect to PV, in cm (with sign)"),
    NTupleVariable("edxy",  lambda x : x.dxyError(), help="#sigma(d_{xy}) with respect to PV, in cm"),
    NTupleVariable("edz", lambda x : x.dzError(), help="#sigma(d_{z}) with respect to PV, in cm"),    

    NTupleVariable("trackerLayers", lambda x : x.hitPattern().trackerLayersWithMeasurement(), int, help="Tracker Layers"),
    NTupleVariable("pixelLayers", lambda x : x.hitPattern().pixelLayersWithMeasurement(), int, help="Pixel Layers"),
    NTupleVariable("trackerHits", lambda x : x.hitPattern().numberOfValidTrackerHits(), int, help="Tracker hits"),
    NTupleVariable("pixelHits", lambda x : x.hitPattern().numberOfValidPixelHits(), int, help="Pixel hits"),
    NTupleVariable("missingInnerPixelHits", lambda x : x.hitPattern().numberOfLostPixelHits(1), int, help="Missing inner pixel hits"),
    NTupleVariable("missingOuterPixelHits", lambda x : x.hitPattern().numberOfLostPixelHits(2), int, help="Missing outer pixel hits"),
    NTupleVariable("missingInnerStripHits", lambda x : x.hitPattern().numberOfLostStripHits(1), int, help="Missing inner strips hits"),
    NTupleVariable("missingOuterStripHits", lambda x : x.hitPattern().numberOfLostStripHits(2), int, help="Missing outer strips hits"),
    NTupleVariable("missingInnerTrackerHits", lambda x : x.hitPattern().numberOfLostTrackerHits(1), int, help="Missing inner tracker hits"),
    NTupleVariable("missingOuterTrackerHits", lambda x : x.hitPattern().numberOfLostTrackerHits(2), int, help="Missing outer tracker hits"),
    NTupleVariable("missingMiddleTrackerHits", lambda x : x.hitPattern().numberOfLostTrackerHits(0), int, help="Missing tracker hits in the middle of the track"),

    NTupleVariable("highPurity", lambda x : x.isHighPurityTrack(), int, help="High purity"),

    NTupleVariable("miniIsoCH",   lambda x : x.miniPFIsolation().chargedHadronIso(), help="Charged hadron mini-isolation"),
    NTupleVariable("miniIsoNH",   lambda x : x.miniPFIsolation().neutralHadronIso(), help="Neutral hadron mini-isolation"),
    NTupleVariable("miniIsoPH",   lambda x : x.miniPFIsolation().photonIso(), help="Photon mini-isolation"),
    NTupleVariable("miniIsoPU",   lambda x : x.miniPFIsolation().puChargedHadronIso(), help="Pileup charged hadron mini-isolation"),
    NTupleVariable("miniRelIso",   lambda x : _isoDBeta(x.miniPFIsolation())/x.pt(), help="mini-isolation (relative, PU-corrected using dBeta)"),
    NTupleVariable("dR03IsoCH",   lambda x : x.pfIsolationDR03().chargedHadronIso(), help="Charged hadron isolation dR=0.3"),
    NTupleVariable("dR03IsoNH",   lambda x : x.pfIsolationDR03().neutralHadronIso(), help="Neutral hadron isolation dR=0.3"),
    NTupleVariable("dR03IsoPH",   lambda x : x.pfIsolationDR03().photonIso(), help="Photon isolation = dR=0.3"),
    NTupleVariable("dR03IsoPU",   lambda x : x.pfIsolationDR03().puChargedHadronIso(), help="Pileup charged hadron isolation dR=0.3"),
    NTupleVariable("relIso03",   lambda x : _isoDBeta(x.pfIsolationDR03())/x.pt(), help="relative isolation dR=0.3 (PU-corrected using dBeta)"),
    NTupleVariable("caloEmEnergy",   lambda x : x.matchedCaloJetEmEnergy(), help="Energy in the ECAL behind the track"),
    NTupleVariable("caloHadEnergy",   lambda x : x.matchedCaloJetHadEnergy(), help="Energy in the HCAL behind the track"),
    NTupleVariable("channelsGoodECAL", lambda x : x.channelsGoodECAL, int, help="Flag set to 1 when the track extrapolates to all good ECAL channels"),
    NTupleVariable("channelsGoodHCAL", lambda x : x.channelsGoodHCAL, int, help="Flag set to 1 when the track extrapolates to all good HCAL channels"),

    NTupleVariable("awayJet_idx", lambda x : x.leadAwayJet.index if x.leadAwayJet else -1, int),
    NTupleVariable("awayJet_pt", lambda x : x.leadAwayJet.pt() if x.leadAwayJet else 0),
    NTupleVariable("awayJet_dr", lambda x : deltaR(x, x.leadAwayJet) if x.leadAwayJet else 0),

    NTupleVariable("awayNJet", lambda x : x.awayJetInfo['num'], int),
    NTupleVariable("awayHTJet", lambda x : x.awayJetInfo['ht']),

    NTupleVariable("awayMu_idx", lambda x : x.leadAwayMu.index if x.leadAwayMu else -1, int),
    NTupleVariable("awayMu_dr", lambda x : deltaR(x, x.leadAwayMu) if x.leadAwayMu else 0),
    NTupleVariable("awayMu_mll", lambda x : (x.leadAwayMu.p4()+x.p4()).M() if x.leadAwayMu else 0),

    NTupleVariable("awayEle_idx", lambda x : x.leadAwayEle.index if x.leadAwayEle else -1, int),
    NTupleVariable("awayEle_dr", lambda x : deltaR(x, x.leadAwayEle) if x.leadAwayEle else 0),
    NTupleVariable("awayEle_mll", lambda x : (x.leadAwayEle.p4()+x.p4()).M() if x.leadAwayEle else 0),

    NTupleVariable("closestMu_idx", lambda x : x.closestMu.index if x.closestMu else -1, int),
    NTupleVariable("closestEle_idx", lambda x : x.closestEle.index if x.closestEle else -1, int),
    NTupleVariable("closestTau_idx", lambda x : x.closestTau.index if x.closestTau else -1, int),

    NTupleVariable("trigLepton_idx", lambda x : x.trigLepton.index if getattr(x, 'trigLepton', None) else -1, int),

    NTupleVariable("myDeDx", lambda x : x.myDeDx),
    NTupleVariable("dedxByLayer0", lambda x : x.dedxByLayer[0]),
    NTupleVariable("dedxByLayer1", lambda x : x.dedxByLayer[1]),
    NTupleVariable("dedxByLayer2", lambda x : x.dedxByLayer[2]),
    NTupleVariable("dedxByLayer3", lambda x : x.dedxByLayer[3]),
    NTupleVariable("dedxByLayer4", lambda x : x.dedxByLayer[4]),
    NTupleVariable("dedxByLayer5", lambda x : x.dedxByLayer[5]),
    NTupleVariable("dedxByLayer6", lambda x : x.dedxByLayer[6]),
    NTupleVariable("dedxByLayer7", lambda x : x.dedxByLayer[7]),
    NTupleVariable("dedxByLayer8", lambda x : x.dedxByLayer[8]),
    NTupleVariable("dedxByLayer9", lambda x : x.dedxByLayer[9]),
    NTupleVariable("dedxByLayer10", lambda x : x.dedxByLayer[10]),
    NTupleVariable("dedxByLayer11", lambda x : x.dedxByLayer[11]),
    NTupleVariable("dedxByLayer12", lambda x : x.dedxByLayer[12]),
    NTupleVariable("dedxByLayer13", lambda x : x.dedxByLayer[13]),
    
    NTupleVariable("subDetIdByLayer0", lambda x : x.subDetIdByLayer[0], int),
    NTupleVariable("subDetIdByLayer1", lambda x : x.subDetIdByLayer[1], int),
    NTupleVariable("subDetIdByLayer2", lambda x : x.subDetIdByLayer[2], int),
    NTupleVariable("subDetIdByLayer3", lambda x : x.subDetIdByLayer[3], int),
    NTupleVariable("subDetIdByLayer4", lambda x : x.subDetIdByLayer[4], int),
    NTupleVariable("subDetIdByLayer5", lambda x : x.subDetIdByLayer[5], int),
    NTupleVariable("subDetIdByLayer6", lambda x : x.subDetIdByLayer[6], int),
    NTupleVariable("subDetIdByLayer7", lambda x : x.subDetIdByLayer[7], int),
    NTupleVariable("subDetIdByLayer8", lambda x : x.subDetIdByLayer[8], int),
    NTupleVariable("subDetIdByLayer9", lambda x : x.subDetIdByLayer[9], int),
    NTupleVariable("subDetIdByLayer10", lambda x : x.subDetIdByLayer[10], int),
    NTupleVariable("subDetIdByLayer11", lambda x : x.subDetIdByLayer[11], int),
    NTupleVariable("subDetIdByLayer12", lambda x : x.subDetIdByLayer[12], int),
    NTupleVariable("subDetIdByLayer13", lambda x : x.subDetIdByLayer[13], int),

    NTupleVariable("sizeXbyLayer0", lambda x : x.sizeXbyLayer[0], int),
    NTupleVariable("sizeXbyLayer1", lambda x : x.sizeXbyLayer[1], int),
    NTupleVariable("sizeXbyLayer2", lambda x : x.sizeXbyLayer[2], int),
    NTupleVariable("sizeXbyLayer3", lambda x : x.sizeXbyLayer[3], int),
    NTupleVariable("sizeXbyLayer4", lambda x : x.sizeXbyLayer[4], int),
    NTupleVariable("sizeXbyLayer5", lambda x : x.sizeXbyLayer[5], int),
    NTupleVariable("sizeXbyLayer6", lambda x : x.sizeXbyLayer[6], int),
    NTupleVariable("sizeXbyLayer7", lambda x : x.sizeXbyLayer[7], int),
    NTupleVariable("sizeXbyLayer8", lambda x : x.sizeXbyLayer[8], int),
    NTupleVariable("sizeXbyLayer9", lambda x : x.sizeXbyLayer[9], int),
    NTupleVariable("sizeXbyLayer10", lambda x : x.sizeXbyLayer[10], int),
    NTupleVariable("sizeXbyLayer11", lambda x : x.sizeXbyLayer[11], int),
    NTupleVariable("sizeXbyLayer12", lambda x : x.sizeXbyLayer[12], int),
    NTupleVariable("sizeXbyLayer13", lambda x : x.sizeXbyLayer[13], int),
                                                                                                         
    NTupleVariable("sizeYbyLayer0", lambda x : x.sizeYbyLayer[0], int),
    NTupleVariable("sizeYbyLayer1", lambda x : x.sizeYbyLayer[1], int),
    NTupleVariable("sizeYbyLayer2", lambda x : x.sizeYbyLayer[2], int),
    NTupleVariable("sizeYbyLayer3", lambda x : x.sizeYbyLayer[3], int),
    NTupleVariable("sizeYbyLayer4", lambda x : x.sizeYbyLayer[4], int),
    NTupleVariable("sizeYbyLayer5", lambda x : x.sizeYbyLayer[5], int),
    NTupleVariable("sizeYbyLayer6", lambda x : x.sizeYbyLayer[6], int),
    NTupleVariable("sizeYbyLayer7", lambda x : x.sizeYbyLayer[7], int),
    NTupleVariable("sizeYbyLayer8", lambda x : x.sizeYbyLayer[8], int),
    NTupleVariable("sizeYbyLayer9", lambda x : x.sizeYbyLayer[9], int),
    NTupleVariable("sizeYbyLayer10", lambda x : x.sizeYbyLayer[10], int),
    NTupleVariable("sizeYbyLayer11", lambda x : x.sizeYbyLayer[11], int),
    NTupleVariable("sizeYbyLayer12", lambda x : x.sizeYbyLayer[12], int),
    NTupleVariable("sizeYbyLayer13", lambda x : x.sizeYbyLayer[13], int),

    # geometry for pixels
    NTupleVariable("layerPixelByLayer0",  lambda x : x.layerPixelByLayer[0] ),
    NTupleVariable("layerPixelByLayer1",  lambda x : x.layerPixelByLayer[1] ),
    NTupleVariable("layerPixelByLayer2",  lambda x : x.layerPixelByLayer[2] ),
    NTupleVariable("layerPixelByLayer3",  lambda x : x.layerPixelByLayer[3] ),
    NTupleVariable("layerPixelByLayer4",  lambda x : x.layerPixelByLayer[4] ),
    NTupleVariable("layerPixelByLayer5",  lambda x : x.layerPixelByLayer[5] ),
    NTupleVariable("layerPixelByLayer6",  lambda x : x.layerPixelByLayer[6] ),
    NTupleVariable("layerPixelByLayer7",  lambda x : x.layerPixelByLayer[7] ),
    NTupleVariable("layerPixelByLayer8",  lambda x : x.layerPixelByLayer[8] ),
    NTupleVariable("layerPixelByLayer9",  lambda x : x.layerPixelByLayer[9] ),
    NTupleVariable("layerPixelByLayer10", lambda x : x.layerPixelByLayer[10]),
    NTupleVariable("layerPixelByLayer11", lambda x : x.layerPixelByLayer[11]),
    NTupleVariable("layerPixelByLayer12", lambda x : x.layerPixelByLayer[12]),
    NTupleVariable("layerPixelByLayer13", lambda x : x.layerPixelByLayer[13]),

    NTupleVariable("ladderPixelByLayer0",  lambda x : x.ladderPixelByLayer[0] ),
    NTupleVariable("ladderPixelByLayer1",  lambda x : x.ladderPixelByLayer[1] ),
    NTupleVariable("ladderPixelByLayer2",  lambda x : x.ladderPixelByLayer[2] ),
    NTupleVariable("ladderPixelByLayer3",  lambda x : x.ladderPixelByLayer[3] ),
    NTupleVariable("ladderPixelByLayer4",  lambda x : x.ladderPixelByLayer[4] ),
    NTupleVariable("ladderPixelByLayer5",  lambda x : x.ladderPixelByLayer[5] ),
    NTupleVariable("ladderPixelByLayer6",  lambda x : x.ladderPixelByLayer[6] ),
    NTupleVariable("ladderPixelByLayer7",  lambda x : x.ladderPixelByLayer[7] ),
    NTupleVariable("ladderPixelByLayer8",  lambda x : x.ladderPixelByLayer[8] ),
    NTupleVariable("ladderPixelByLayer9",  lambda x : x.ladderPixelByLayer[9] ),
    NTupleVariable("ladderPixelByLayer10", lambda x : x.ladderPixelByLayer[10]),
    NTupleVariable("ladderPixelByLayer11", lambda x : x.ladderPixelByLayer[11]),
    NTupleVariable("ladderPixelByLayer12", lambda x : x.ladderPixelByLayer[12]),
    NTupleVariable("ladderPixelByLayer13", lambda x : x.ladderPixelByLayer[13]),

    NTupleVariable("modulePixelByLayer0",  lambda x : x.modulePixelByLayer[0] ),
    NTupleVariable("modulePixelByLayer1",  lambda x : x.modulePixelByLayer[1] ),
    NTupleVariable("modulePixelByLayer2",  lambda x : x.modulePixelByLayer[2] ),
    NTupleVariable("modulePixelByLayer3",  lambda x : x.modulePixelByLayer[3] ),
    NTupleVariable("modulePixelByLayer4",  lambda x : x.modulePixelByLayer[4] ),
    NTupleVariable("modulePixelByLayer5",  lambda x : x.modulePixelByLayer[5] ),
    NTupleVariable("modulePixelByLayer6",  lambda x : x.modulePixelByLayer[6] ),
    NTupleVariable("modulePixelByLayer7",  lambda x : x.modulePixelByLayer[7] ),
    NTupleVariable("modulePixelByLayer8",  lambda x : x.modulePixelByLayer[8] ),
    NTupleVariable("modulePixelByLayer9",  lambda x : x.modulePixelByLayer[9] ),
    NTupleVariable("modulePixelByLayer10", lambda x : x.modulePixelByLayer[10]),
    NTupleVariable("modulePixelByLayer11", lambda x : x.modulePixelByLayer[11]),
    NTupleVariable("modulePixelByLayer12", lambda x : x.modulePixelByLayer[12]),
    NTupleVariable("modulePixelByLayer13", lambda x : x.modulePixelByLayer[13]),

    NTupleVariable("sidePixelByLayer0",  lambda x : x.sidePixelByLayer[0] ),
    NTupleVariable("sidePixelByLayer1",  lambda x : x.sidePixelByLayer[1] ),
    NTupleVariable("sidePixelByLayer2",  lambda x : x.sidePixelByLayer[2] ),
    NTupleVariable("sidePixelByLayer3",  lambda x : x.sidePixelByLayer[3] ),
    NTupleVariable("sidePixelByLayer4",  lambda x : x.sidePixelByLayer[4] ),
    NTupleVariable("sidePixelByLayer5",  lambda x : x.sidePixelByLayer[5] ),
    NTupleVariable("sidePixelByLayer6",  lambda x : x.sidePixelByLayer[6] ),
    NTupleVariable("sidePixelByLayer7",  lambda x : x.sidePixelByLayer[7] ),
    NTupleVariable("sidePixelByLayer8",  lambda x : x.sidePixelByLayer[8] ),
    NTupleVariable("sidePixelByLayer9",  lambda x : x.sidePixelByLayer[9] ),
    NTupleVariable("sidePixelByLayer10", lambda x : x.sidePixelByLayer[10]),
    NTupleVariable("sidePixelByLayer11", lambda x : x.sidePixelByLayer[11]),
    NTupleVariable("sidePixelByLayer12", lambda x : x.sidePixelByLayer[12]),
    NTupleVariable("sidePixelByLayer13", lambda x : x.sidePixelByLayer[13]),

    NTupleVariable("diskPixelByLayer0",  lambda x : x.diskPixelByLayer[0] ),
    NTupleVariable("diskPixelByLayer1",  lambda x : x.diskPixelByLayer[1] ),
    NTupleVariable("diskPixelByLayer2",  lambda x : x.diskPixelByLayer[2] ),
    NTupleVariable("diskPixelByLayer3",  lambda x : x.diskPixelByLayer[3] ),
    NTupleVariable("diskPixelByLayer4",  lambda x : x.diskPixelByLayer[4] ),
    NTupleVariable("diskPixelByLayer5",  lambda x : x.diskPixelByLayer[5] ),
    NTupleVariable("diskPixelByLayer6",  lambda x : x.diskPixelByLayer[6] ),
    NTupleVariable("diskPixelByLayer7",  lambda x : x.diskPixelByLayer[7] ),
    NTupleVariable("diskPixelByLayer8",  lambda x : x.diskPixelByLayer[8] ),
    NTupleVariable("diskPixelByLayer9",  lambda x : x.diskPixelByLayer[9] ),
    NTupleVariable("diskPixelByLayer10", lambda x : x.diskPixelByLayer[10]),
    NTupleVariable("diskPixelByLayer11", lambda x : x.diskPixelByLayer[11]),
    NTupleVariable("diskPixelByLayer12", lambda x : x.diskPixelByLayer[12]),
    NTupleVariable("diskPixelByLayer13", lambda x : x.diskPixelByLayer[13]),

    NTupleVariable("bladePixelByLayer0",  lambda x : x.bladePixelByLayer[0] ),
    NTupleVariable("bladePixelByLayer1",  lambda x : x.bladePixelByLayer[1] ),
    NTupleVariable("bladePixelByLayer2",  lambda x : x.bladePixelByLayer[2] ),
    NTupleVariable("bladePixelByLayer3",  lambda x : x.bladePixelByLayer[3] ),
    NTupleVariable("bladePixelByLayer4",  lambda x : x.bladePixelByLayer[4] ),
    NTupleVariable("bladePixelByLayer5",  lambda x : x.bladePixelByLayer[5] ),
    NTupleVariable("bladePixelByLayer6",  lambda x : x.bladePixelByLayer[6] ),
    NTupleVariable("bladePixelByLayer7",  lambda x : x.bladePixelByLayer[7] ),
    NTupleVariable("bladePixelByLayer8",  lambda x : x.bladePixelByLayer[8] ),
    NTupleVariable("bladePixelByLayer9",  lambda x : x.bladePixelByLayer[9] ),
    NTupleVariable("bladePixelByLayer10", lambda x : x.bladePixelByLayer[10]),
    NTupleVariable("bladePixelByLayer11", lambda x : x.bladePixelByLayer[11]),
    NTupleVariable("bladePixelByLayer12", lambda x : x.bladePixelByLayer[12]),
    NTupleVariable("bladePixelByLayer13", lambda x : x.bladePixelByLayer[13]),

    NTupleVariable("panelPixelByLayer0",  lambda x : x.panelPixelByLayer[0] ),
    NTupleVariable("panelPixelByLayer1",  lambda x : x.panelPixelByLayer[1] ),
    NTupleVariable("panelPixelByLayer2",  lambda x : x.panelPixelByLayer[2] ),
    NTupleVariable("panelPixelByLayer3",  lambda x : x.panelPixelByLayer[3] ),
    NTupleVariable("panelPixelByLayer4",  lambda x : x.panelPixelByLayer[4] ),
    NTupleVariable("panelPixelByLayer5",  lambda x : x.panelPixelByLayer[5] ),
    NTupleVariable("panelPixelByLayer6",  lambda x : x.panelPixelByLayer[6] ),
    NTupleVariable("panelPixelByLayer7",  lambda x : x.panelPixelByLayer[7] ),
    NTupleVariable("panelPixelByLayer8",  lambda x : x.panelPixelByLayer[8] ),
    NTupleVariable("panelPixelByLayer9",  lambda x : x.panelPixelByLayer[9] ),
    NTupleVariable("panelPixelByLayer10", lambda x : x.panelPixelByLayer[10]),
    NTupleVariable("panelPixelByLayer11", lambda x : x.panelPixelByLayer[11]),
    NTupleVariable("panelPixelByLayer12", lambda x : x.panelPixelByLayer[12]),
    NTupleVariable("panelPixelByLayer13", lambda x : x.panelPixelByLayer[13]),

                                                                                                       
    NTupleVariable("mcMatch", lambda x : x.mcMatch.index if x.mcMatch else -1, int, mcOnly=True),
    NTupleVariable("mcMatchAnyId", lambda x : x.mcMatchAny.pdgId()*(1+99*x.mcMatchAny.isDirectPromptTauDecayProductFinalState()) if x.mcMatchAny else 0, int, mcOnly=True, help="MC pdgId of the matched gen lepton, tau, photon or chargino (for leptons from tau, it's pdgId*100)"),
    NTupleVariable("mcMatchAnyPt", lambda x : x.mcMatchAny.pt() if x.mcMatchAny else 0, int, mcOnly=True, help="MC pt of the matched gen lepton, tau, photon or chargino"),
]
)

##------------------------------------------  
## genCharginoType
##------------------------------------------  
genCharginoType  = NTupleObjectType("genCharginoType", baseObjectTypes = [ genParticleWithMotherId ], variables = [
    NTupleVariable("beta", lambda x : x.p()/x.energy()),
    NTupleVariable("decayR", lambda x : x.decayPoint.R()),
    NTupleVariable("decayZ", lambda x : x.decayPoint.Z()),
])

## Tree Producer
treeProducer = cfg.Analyzer(
    AutoFillTreeProducer, name='treeProducerXtracks',
    vectorTree = True, saveTLorentzVectors = False,  defaultFloatType = 'F', PDFWeights = [],
    globalVariables = [
        NTupleVariable("rho",  lambda ev: ev.rho, float, help="kt6PFJets rho"),
        NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"),

        NTupleVariable("nJet30", lambda ev: sum([j.pt() > 30 for j in ev.cleanJets]), int, help="Number of jets with pt > 30, |eta|<2.4"),
        NTupleVariable("nJet30a", lambda ev: sum([j.pt() > 30 for j in ev.cleanJetsAll]), int, help="Number of jets with pt > 30, |eta|<4.7"),

        ## ------- lheHT, needed for merging HT binned samples
        NTupleVariable("lheHT", lambda ev : getattr(ev,"lheHT",-999), mcOnly=True, help="H_{T} computed from quarks and gluons in Heppy LHEAnalyzer"),
        NTupleVariable("lheHTIncoming", lambda ev : getattr(ev,"lheHTIncoming",-999), mcOnly=True, help="H_{T} computed from quarks and gluons in Heppy LHEAnalyzer (only LHE status<0 as mothers)"),
        NTupleVariable("lheVpt", lambda ev : getattr(ev,"lheV_pt",-999), mcOnly=True, help="p_{T} of the V boson at LHE level"),
        ],
    globalObjects = {
        "met"             : NTupleObject("met", metTypeXtracks, help="PF E_{T}^{miss}"),
        "metNoMu"         : NTupleObject("metNoMu", metTypeXtracksBasic, help="PF E_{T}^{miss} without muons"),
        "met_jecUp"       : NTupleObject("met_jecUp", metTypeXtracks, help="PF E_{T}^{miss}, after type 1 corrections (JEC plus 1sigma)"),
        "metNoMu_jecUp"   : NTupleObject("metNoMu_jecUp", metTypeXtracksBasic, help="PF E_{T}^{miss} without muons, after type 1 corrections (JEC plus 1sigma)"),
        "met_jecDown"     : NTupleObject("met_jecDown", metTypeXtracks, help="PF E_{T}^{miss}, after type 1 corrections (JEC minus 1sigma)"),
        "metNoMu_jecDown" : NTupleObject("metNoMu_jecUp", metTypeXtracksBasic, help="PF E_{T}^{miss} without muons, after type 1 corrections (JEC minus 1sigma)"),
        },
    collections = {
        ##--------------------------------------------------
        "selectedTaus"    : NTupleCollection("TauGood",  tauTypeXtracks, 8, help="Taus after the preselection"),
        "selectedLeptons" : NTupleCollection("LepGood",  leptonTypeXtracks, 8, help="Leptons after the preselection"),
        ##------------------------------------------------
        "cleanJets"       : NTupleCollection("Jet",     jetTypeXtracks, 15, help="Cental jets after full selection and cleaning, sorted by pt"),
        "cleanJetsFwd"    : NTupleCollection("JetFwd",  jetTypeXtracksFwd,  6, help="Forward jets after full selection and cleaning, sorted by pt"),
        ##------------------------------------------------
        "isoTracks"       : NTupleCollection("IsoTrack",  isoTrackTypeDeDx, 4, help="Isolated tracks"),
        ##------------------------------------------------
        "genCharginos"    : NTupleCollection("GenChargino",  genCharginoType, 4, mcOnly=True, help="Gen chargino"),
        },
    )

