from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.Electron import Electron

import ROOT

from PhysicsTools.HeppyCore.utils.deltar import deltaR, matchObjectCollection3

def leading(objs):
    if len(objs) == 0: return None
    return max(objs, key = lambda obj : obj.pt())
def closest(center, objs):
    if len(objs) == 0: return None
    return min(objs, key = lambda obj : deltaR(center.eta(), center.phi(), obj.eta(), obj.phi()))
def nearby(center, objs, dr):
    return [ o for o in objs if deltaR(center.eta(), center.phi(), o.eta(), o.phi()) < dr ]
def cleaned(center, objs, dr):
    return [ o for o in objs if deltaR(center.eta(), center.phi(), o.eta(), o.phi()) >= dr ]
def jetStats(objs):
    return dict( ht = sum((o.pt() for o in objs), 0),
                 num = len(objs) )

class isoTrackDeDxAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(isoTrackDeDxAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)

    def declareHandles(self):
        super(isoTrackDeDxAnalyzer, self).declareHandles()
        self.handles['dedx'] = AutoHandle("isolatedTracks","edm::Association<reco::DeDxHitInfoCollection>")            

    def beginLoop(self, setup):
        super(isoTrackDeDxAnalyzer,self).beginLoop(setup)
        self.counters.addCounter('events')
        count = self.counters.counter('events')
        count.register('all events')
        count.register('selected iso track')
        count.register('accepted events')

    def process(self, event):
        self.readCollections( event.input )
        self.counters.counter('events').inc('all events')

        muons = []
        electrons = []
        for i,l in enumerate(event.selectedLeptons):
            l.index = i
            if abs(l.pdgId()) == 13: muons.append(l)
            if abs(l.pdgId()) == 11: electrons.append(l)

        for i,l in enumerate(event.selectedTaus):
            l.index = i

        for i,l in enumerate(event.cleanJets):
            l.index = i

        getDeDxRef = self.handles['dedx'].product().get
        if self.cfg_ana.doDeDx == "94XMiniAODv1-Hack":
            fixed = self.handles['dedx'].product().fixOffsets(-1)
            getDeDxRef = fixed.get
        
        event.isoTracks = []
        
        pixelChargeToEnergyCoefficient = 3.61e-6
        stripChargeToEnergyCoefficient = 3.61e-6 * 265
        
        for t in event.preselIsoTracks:
            # add more variables
            t.leadAwayJet = leading(cleaned(t,  event.cleanJets, 0.4))
            t.awayJetInfo = jetStats(cleaned(t, event.cleanJets, 0.4))
            t.leadAwayMu  = leading(cleaned(t, muons, 0.4))
            t.leadAwayEle = leading(cleaned(t, electrons, 0.4))
            t.leadAwayTau = leading(cleaned(t, event.selectedTaus, 0.4))
            t.closestMu   = closest(t, nearby(t, muons, 0.4))
            t.closestEle  = closest(t, nearby(t, electrons, 0.4))
            t.closestTau  = closest(t, nearby(t, event.selectedTaus, 0.4))

            dedxArray = []
            subDetIdArray = []
            sizeXarray = []
            sizeYarray = []
            
            # geometry for pixels
            layerPixelArray  = []
            ladderPixelArray = []
            modulePixelArray = []
            sidePixelArray   = []
            diskPixelArray   = []
            bladePixelArray  = []
            panelPixelArray  = []

            
            for i in xrange(100):
              dedxArray.append(0)
              subDetIdArray.append(0)
              sizeXarray.append(0)
              sizeYarray.append(0)

              # geometry for pixels              
              layerPixelArray  .append(0)
              ladderPixelArray .append(0)
              modulePixelArray .append(0)
              sidePixelArray   .append(0)
              diskPixelArray   .append(0)
              bladePixelArray  .append(0)
              panelPixelArray  .append(0)
              

            # get dedx
            if self.cfg_ana.doDeDx:
                ref = getDeDxRef(t.index)
                if ref.isNull():
                    print "ERROR: no dE/dx for track of pt %.2f, eta %.2f" % (t.pt(),t.eta())
                    continue
                dedx = ref.get(); 
                nhits = dedx.size()
                
                # this below is just dummy to give you a template
                mysum = 0
                
                for ih in xrange(nhits):
                    pixelCluster = dedx.pixelCluster(ih)
                    stripCluster = dedx.stripCluster(ih)
                    if pixelCluster:
                      dedxArray[ih] = pixelCluster.charge()/dedx.pathlength(ih)
                      # convert number of electrons to MeV
                      dedxArray[ih] *= pixelChargeToEnergyCoefficient
                      
                      sizeXarray[ih] = pixelCluster.sizeX()
                      sizeYarray[ih] = pixelCluster.sizeY()
                      
                      mysum += pixelCluster.charge()

                      #
                      #  https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/DeDxHitInfo.h
                      #
                      #  dEdX object DeDxHitInfoCollection
                      #
                      #
                      #  https://github.com/cms-sw/cmssw/blob/master/DataFormats/SiPixelDetId/interface/PXBDetId.h
                      #  https://github.com/cms-sw/cmssw/blob/master/DataFormats/SiPixelDetId/interface/PXFDetId.h
                      #
                      if dedx.detId(ih).subdetId() == 1 :
                        # barrel
                        layerPixelArray  [ih] = ROOT.PXBDetId( dedx.detId(ih).rawId() ).layer()
                        ladderPixelArray [ih] = ROOT.PXBDetId( dedx.detId(ih).rawId() ).ladder()
                        modulePixelArray [ih] = ROOT.PXBDetId( dedx.detId(ih).rawId() ).module()
                        sidePixelArray   [ih] = -99
                        diskPixelArray   [ih] = -99
                        bladePixelArray  [ih] = -99
                        panelPixelArray  [ih] = -99

                        #print " Barrel[", dedx.detId(ih).subdetId(), "]:: ", ROOT.PXBDetId( dedx.detId(ih).rawId() ).layer(), "   ",  ROOT.PXBDetId( dedx.detId(ih).rawId() ).ladder(), "   ",  ROOT.PXBDetId( dedx.detId(ih).rawId() ).module()                        
                      else :
                        # endcap
                        layerPixelArray  [ih] = -99
                        ladderPixelArray [ih] = -99
                        modulePixelArray [ih] = ROOT.PXFDetId( dedx.detId(ih).rawId() ).module()
                        sidePixelArray   [ih] = ROOT.PXFDetId( dedx.detId(ih).rawId() ).side()
                        diskPixelArray   [ih] = ROOT.PXFDetId( dedx.detId(ih).rawId() ).disk()
                        bladePixelArray  [ih] = ROOT.PXFDetId( dedx.detId(ih).rawId() ).blade()
                        panelPixelArray  [ih] = ROOT.PXFDetId( dedx.detId(ih).rawId() ).panel()

                        #print " Endcap[", dedx.detId(ih).subdetId(), "]:: ", ROOT.PXFDetId( dedx.detId(ih).rawId() ).side(), "   ",  ROOT.PXFDetId( dedx.detId(ih).rawId() ).disk(), "   ",  ROOT.PXFDetId( dedx.detId(ih).rawId() ).blade(), "   ",  ROOT.PXFDetId( dedx.detId(ih).rawId() ).panel(), "   ",  ROOT.PXFDetId( dedx.detId(ih).rawId() ).module()
                        
                        
                    if stripCluster:
                      dedxArray[ih] = stripCluster.charge()/dedx.pathlength(ih)
                      # convert number of electrons to MeV
                      dedxArray[ih] *= stripChargeToEnergyCoefficient
                    subDetIdArray[ih] = dedx.detId(ih).subdetId()

                    
                t.myDeDx = mysum
            else:
                t.myDeDx = 0


            t.dedxByLayer = dedxArray
            t.subDetIdByLayer = subDetIdArray
            t.sizeXbyLayer = sizeXarray
            t.sizeYbyLayer = sizeYarray

            t.layerPixelByLayer   = layerPixelArray  
            t.ladderPixelByLayer  = ladderPixelArray 
            t.modulePixelByLayer  = modulePixelArray 
            t.sidePixelByLayer    = sidePixelArray   
            t.diskPixelByLayer    = diskPixelArray   
            t.bladePixelByLayer   = bladePixelArray  
            t.panelPixelByLayer   = panelPixelArray  

            # add a flag for bad ECAL channels in the way of the track
            t.channelsGoodECAL = 1
            for ie in t.crossedEcalStatus():
                if ie != 0: t.channelsGoodECAL = 0

            # add a flag for bad HCAL channels in the way of the track
            t.channelsGoodHCAL = 1
            for ih in t.crossedHcalStatus():
                if (ih & (1<<5)) != 0: t.channelsGoodHCAL = 0


            # add to the list
            event.isoTracks.append(t)

        if len(event.isoTracks) == 0: 
            return False
        self.counters.counter('events').inc('selected iso track')

        if self.cfg_comp.isMC: # do MC matching
            event.genCharginos = [ g for g in event.generatorSummary if abs(g.pdgId()) in (1000024,1000037)]
            for i,g in enumerate(event.genCharginos):
                g.index = i
                g.decayPoint = g.vertex()
                if g.numberOfDaughters() > 1:
                    for i in xrange(g.numberOfDaughters()):
                        dau = g.daughter(i)
                        if dau: 
                            g.decayPoint = dau.vertex()
                            break
                #print "GEN pdgId %+8d pt %7.1f eta %+6.2f phi %+6.2f status %3d daughters %d decay R %5.1f, z %+5.1f" % (g.pdgId(), g.pt(), g.eta(), g.phi(), g.status(), g.numberOfDaughters(), g.decayPoint.R(), g.decayPoint.Z())
            match = matchObjectCollection3(event.isoTracks, event.genCharginos, deltaRMax = 0.1)
            for t in event.isoTracks:
                t.mcMatch = match[t]
                #print "REC charge %+7d pt %7.1f eta %+6.2f phi %+6.2f" % (t.charge(), t.pt(), t.eta(), t.phi())
            #    for t in event.isoTracks:
            #        if deltaR(g,t) < 0.3:
            #            print " -> charge %+7d pt %7.1f eta %+6.2f phi %+6.2f   dr %.4f" % (t.charge(), t.pt(), t.eta(), t.phi(), deltaR(g,t))
            #print "\n"

            ## Now we add a generic match to charginos + prompt leptons + prompt taus + prompt photons
            anyMatchable = event.genCharginos[:]
            anyMatchable += [ x for x in event.genParticles if abs(x.pdgId()) in (11,13) and (x.isPromptFinalState() or x.isDirectPromptTauDecayProductFinalState()) ]
            anyMatchable += [ x for x in event.genParticles if abs(x.pdgId()) == 15 and x.isPromptDecayed() and x.pt() > 10 ]
            anyMatchable += [ x for x in event.genParticles if x.pdgId() == 22 and x.isPromptFinalState() and x.pt() > 20 ]
            matchAny = matchObjectCollection3(event.isoTracks, anyMatchable, deltaRMax = 0.2)
            for t in event.isoTracks:
                t.mcMatchAny = matchAny[t]


        # do any more event selection
        self.counters.counter('events').inc('accepted events')
        return True
