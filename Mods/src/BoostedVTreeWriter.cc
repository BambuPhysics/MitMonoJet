#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitMonoJet/Mods/interface/BoostedVTreeWriter.h"

using namespace mithep;

ClassImp(mithep::BoostedVTreeWriter)

//--------------------------------------------------------------------------------------------------
BoostedVTreeWriter::BoostedVTreeWriter(const char *name, const char *title) : 
  BaseMod                (name,title),
  fIsData                (kTRUE),
  fMcPartsName           (Names::gkMCPartBrn),
  fMcParts               (0),
  fTriggerObjsName       ("HltObjsMonoJet"),
  fTrigObjs              (0),
  fJetsName              (Names::gkPFJetBrn),
  fJetsFromBranch        (kTRUE),
  fJets                  (0),
  fPfCandidatesName      (Names::gkPFCandidatesBrn),
  fPfCandidatesFromBranch(kTRUE),
  fPfCandidates          (0),
  fPhotonsName           (Names::gkPhotonBrn),
  fPhotonsFromBranch     (kTRUE),
  fPhotons               (0),
  fPfTausName            (Names::gkPFTauBrn),
  fPfTausFromBranch      (kTRUE),
  fPfTaus                (0),
  fLeptonsName           (ModNames::gkMergedLeptonsName),
  fPfNoPileUpName        ("pfnopileupcands"),
  fPfPileUpName          ("pfpileupcands"),
  fQgTaggerCHS           (kTRUE),
  fPileUpDenName         (Names::gkPileupEnergyDensityBrn),
  fPileUpDen             (0),
  fVertexesName          (ModNames::gkGoodVertexesName),
  fJetTriggerObjs        (0),
  fConeSize              (0.6),
  fPrune                 (1.),
  fNAnalyzed             (0),
  fOutputName            ("BoostedVNtuple.root"),
  fOutputFile            (0)
{
  // Constructor.

  // WARNING, defining the object here invalidates the call of the setter for the CHS flag
  fQgTagger = new QGTagger(fQgTaggerCHS);
}

BoostedVTreeWriter::~BoostedVTreeWriter()
{
  // Destructor
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void BoostedVTreeWriter::Process()
{
  // Load the branches we want to work with
  LoadEventObject(fJetsName,fJets,fJetsFromBranch);
  LoadEventObject(fPfCandidatesName,fPfCandidates,fPfCandidatesFromBranch);
  LoadEventObject(fPhotonsName,fPhotons,fPhotonsFromBranch);
  LoadEventObject(fPfTausName,fPfTaus,fPfTausFromBranch);
  LoadEventObject(fPfTausName,fPfTaus,fPfTausFromBranch);
  LoadEventObject(fPileUpDenName,fPileUpDen,kTRUE);
  const ParticleOArr *leptons = static_cast<const ParticleOArr*>(GetObject<ParticleOArr>(ModNames::gkMergedLeptonsName));

  // Careful the are not booked and just local (needed only for isolation)
  const PFCandidateCol *lPFNoPileUpCands = GetObject<PFCandidateCol>(fPfNoPileUpName);    
  const PFCandidateCol *lPFPileUpCands = GetObject<PFCandidateCol>(fPfPileUpName);

  // Also these are not booked and just local (needed only for QG discriminant)
  const VertexCol *lVertexes = GetObject<VertexOArr>(fVertexesName);

  // Extract the jet trigger objects from all trigger objects
  GetJetTriggerObjs();

  // Initializes all variables
  fMitGPTree.InitVariables();

  // After tree is initialized we can start with MC if applicable
  if (! fIsData)
    ProcessMc();
 
  // Keep track of events analyzed
  fNAnalyzed++;

  fMitGPTree.run_ = GetEventHeader()->RunNum();
  fMitGPTree.event_ = GetEventHeader()->EvtNum();
  fMitGPTree.lumi_ = GetEventHeader()->LumiSec();
  
  // Tell to the QG tagger what is the PU energy density
  fQgTagger->SetRhoIso(fPileUpDen->At(0)->RhoRandomLowEta());
  float jet1QgTag = 0;
  float jet2QgTag = 0;
  
  // Loop over jets and perform Nsubjettiness analysis (for now just stick with the first jet)
  std::vector<fastjet::PseudoJet> lFjParts;
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {      
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(i));
    if (! jet) {
      printf(" BoostedVTreeWriter::Process() - ERROR - jets provided are not PFJets.");
      break;
    }
    
    // Figure out whether the jet is matched with one of the triggers
    fTrigObjs = GetHLTObjects(fTriggerObjsName);
    if (! fTrigObjs)
      printf(" BoostedVTreeWriter::Process() - ERROR - TriggerObjectCol not found\n");
    else {
      // loop through the stored trigger objects and find corresponding trigger name
      for (UInt_t j=0; j<fTrigObjs->GetEntries();++j) {
    	const TriggerObject *to = fTrigObjs->At(j);
    	TString trName = to->TrigName();
    	// default MonoJet
    	if (trName.Contains("MonoCentralPFJet80_PFMETnoMu"))
    	  fMitGPTree.trigger_ |= 1 << 0;
	if (trName.Contains("HLT_MET120_HBHENoiseCleaned_v"))
    	  fMitGPTree.trigger_ |= 1 << 1;
      }
    }

    // Compute the QG discrimination varaible for this jet and assign it
    fQgTagger->CalculateVariables(jet, lVertexes);
    if       (i == 0)
      jet1QgTag = fQgTagger->QGValue();
    else if (i == 1)
      jet2QgTag = fQgTagger->QGValue();    
    
    // Push all particle flow candidates into fastjet particle collection
    for (UInt_t j=0; j<jet->NPFCands(); ++j) {      
      const PFCandidate *pfCand = jet->PFCand(j);
      lFjParts.push_back(fastjet::PseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E()));
      lFjParts.back().set_user_index(j);
    }	
    break; // this is for now just considering the first jet (consider all in the future)
  }

  // ---- Fastjet is ready ----

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea *lClustering =
    new fastjet::ClusterSequenceArea(lFjParts,*fCAJetDef,*fAreaDefinition);

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(lClustering->inclusive_jets(0.0));
 
  // Access the 2 hardest jets as applicable
  fastjet::PseudoJet lJet1, lJet2;
  if (fPrune) {
    if (lOutJets.size() > 0)
      lJet1 = (*fPruner)(lOutJets[0]);
    if (lOutJets.size() > 1)
      lJet2 = (*fPruner)(lOutJets[1]);
  }
  else {
    if (lOutJets.size() > 0)
      lJet1 = lOutJets[0];
    if (lOutJets.size() > 1)
      lJet2 = lOutJets[1];
  }

  // ---- Fastjet is already done ----

  // Things to monitor before cuts
  fMitGPTree.nParts_ = fPfCandidates->GetEntries();
  fMitGPTree.numJets_ = lOutJets.size();

  // Skip event if no jets are found
  if (lOutJets.size() == 0) {
    if (lClustering)
      delete lClustering;
    return;
  }
    
  // Remove events where hardest jet is lower than 100 GeV
  if (lJet1.pt() < 100) {
    if (lClustering)
      delete lClustering;
    return;
  }
  
  // Fill up the tree
  fMitGPTree.jet1_.SetPxPyPzE( lJet1.px(), lJet1.py(), lJet1.pz(), lJet1.e() );
  fMitGPTree.jet1NParts_ = lJet1.constituents().size();
  fMitGPTree.jet1Pt_ = lJet1.pt();
  fMitGPTree.jet1Eta_ = lJet1.eta();
  fMitGPTree.jet1Phi_ = lJet1.phi();
  fMitGPTree.jet1M_ = lJet1.m();
  fMitGPTree.jet1Tau1_ = GetTau(lJet1,1,1);
  fMitGPTree.jet1Tau2_ = GetTau(lJet1,2,1);
  fMitGPTree.jet1Tau3_ = GetTau(lJet1,3,1);
  fMitGPTree.jet1MinTrigDr_ = MinTriggerDeltaR(fMitGPTree.jet1_);
  fMitGPTree.jet1QgTag_ = jet1QgTag;

  // Only if there is a second jet
  if (lOutJets.size() > 1) {
    fMitGPTree.jet2_.SetPxPyPzE(lJet2.px(),lJet2.py(),lJet2.pz(),lJet2.e());
    fMitGPTree.jet2NParts_ = lJet2.constituents().size();
    fMitGPTree.jet2Pt_ = lJet2.pt();
    fMitGPTree.jet2Eta_ = lJet2.eta();
    fMitGPTree.jet2Phi_ = lJet2.phi();
    fMitGPTree.jet2M_ = lJet2.m();
    fMitGPTree.jet2Tau1_ = GetTau(lJet2,1,1);
    fMitGPTree.jet2Tau2_ = GetTau(lJet2,2,1);
    fMitGPTree.jet2Tau3_ = GetTau(lJet2,3,1);
    fMitGPTree.jet2MinTrigDr_ = MinTriggerDeltaR(fMitGPTree.jet2_);
    fMitGPTree.jet2QgTag_ = jet2QgTag;
  }

  // LEPTONS (MUONS/ELECTRONS)

  fMitGPTree.nLep_ = leptons->GetEntries();
  if (leptons->GetEntries() >= 1) {         // loop over all leptons
    const Particle *lep = leptons->At(0);
    
    fMitGPTree.lep1_ = lep->Mom();
    if      (lep->ObjType() == kMuon) {
      fMitGPTree.lepId1_ = 13;
      const Muon* mu = dynamic_cast<const Muon*>(lep);
      fMitGPTree.lep1IsTightMuon_ = IsTightMuon(mu);
      fMitGPTree.lep1PtErr_ = mu->BestTrk()->PtErr()/mu->BestTrk()->Pt();
      double totalIso =
	IsolationTools::BetaMwithPUCorrection(lPFNoPileUpCands,lPFPileUpCands,mu,0.4);
      fMitGPTree.lep1IsIsolated_ = totalIso < (mu->Pt()*0.2);
    }
    else if (lep->ObjType() == kElectron)
      fMitGPTree.lepId1_ = 11;
    else
      assert(0);  // cannot happen: leptons are only muons and electrons

    if (lep->Charge() < 0)
      fMitGPTree.lepId1_ = -1 * fMitGPTree.lepId1_;
  }

  if (leptons->GetEntries() >= 2) {
    const Particle *lep = leptons->At(1);

    fMitGPTree.lep2_ = lep->Mom();
    if     (lep->ObjType() == kMuon) {
      fMitGPTree.lepId2_ = 13;
      const Muon* mu = dynamic_cast<const Muon*>(lep);
      fMitGPTree.lep2IsTightMuon_ = IsTightMuon(mu);
      fMitGPTree.lep2PtErr_ = mu->BestTrk()->PtErr()/mu->BestTrk()->Pt();
      double totalIso =
	IsolationTools::BetaMwithPUCorrection(lPFNoPileUpCands,lPFPileUpCands,mu,0.4);
      fMitGPTree.lep2IsIsolated_ = totalIso < (mu->Pt()*0.2);
    }
    else if (lep->ObjType() == kElectron)
      fMitGPTree.lepId2_ = 11;
    else
      assert(0);

    if (lep->Charge() < 0)
      fMitGPTree.lepId2_ = -1 * fMitGPTree.lepId2_;
  }

  if (leptons->GetEntries() >= 3) {
    const Particle *lep = leptons->At(2);

    fMitGPTree.lep3_ = lep->Mom();
    if      (lep->ObjType() == kMuon) {
      fMitGPTree.lepId3_ = 13;
      fMitGPTree.lep3IsTightMuon_ = IsTightMuon(dynamic_cast<const Muon*>(leptons->At(2)));
    }
    else if (lep->ObjType() == kElectron)
      fMitGPTree.lepId3_ = 11;
    else
      assert(0);

    if (lep->Charge() < 0)
      fMitGPTree.lepId3_ = -1 * fMitGPTree.lepId3_;

    // If the event contains more than 2 we have no further assumption ( or should we? WZ ;-) )

  }

  // PHOTON(S)

  fMitGPTree.nPhotons_ = fPhotons->GetEntries();
  if (fPhotons->GetEntries() >= 1) {
    const Photon *photon = fPhotons->At(0);
    fMitGPTree.pho1_ = photon->Mom();
  }

  // TAUS

  fMitGPTree.nTaus_ = fPfTaus->GetEntries();
  if (fPfTaus->GetEntries() >= 1) {
    const PFTau *tau = fPfTaus->At(0);
    fMitGPTree.tau1_ = tau->Mom();
    if (fPfTaus->GetEntries() >= 2) {
      tau = fPfTaus->At(1);
      fMitGPTree.tau2_ = tau->Mom();
    }
  }
 
  // Called once when all is done
  fMitGPTree.tree_->Fill();

  // Always cleanup
  if (lClustering)
    delete lClustering;
}

//--------------------------------------------------------------------------------------------------
void BoostedVTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // particle flow collection branch.

  if (! fIsData) {
    ReqEventObject(fMcPartsName,fMcParts,kTRUE);
  }
  ReqEventObject(fJetsName,fJets,fJetsFromBranch);
  ReqEventObject(fPfCandidatesName,fPfCandidates,fPfCandidatesFromBranch);
  ReqEventObject(fPhotonsName,fPhotons,fPhotonsFromBranch);
  ReqEventObject(fPfTausName,fPfTaus,fPfTausFromBranch);
  ReqEventObject(fPileUpDenName,fPileUpDen,kTRUE);

  // Default pruning parameters
  fPruner          = new fastjet::Pruner(fastjet::cambridge_algorithm,0.1,0.5);      // CMS Default
  
  // // CA constructor (fConeSize = 0.8 for CA8) // CMS default
  // fCAJetDef       = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fConeSize);

  // CA constructor (fConeSize = 0.6 for antiKt) - reproducing paper 1: http://arxiv.org/abs/1011.2268
  fCAJetDef       = new fastjet::JetDefinition(fastjet::antikt_algorithm, fConeSize);
  
  // Initialize area caculation (done with ghost particles)
  int    activeAreaRepeats = 1;
  double ghostArea         = 0.01;
  double ghostEtaMax       = 7.0;
  fActiveArea      = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition  = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);

  // Create Ntuple Tree

  fOutputFile = TFile::Open(fOutputName,"RECREATE");
  fMitGPTree.CreateTree();
  fMitGPTree.tree_->SetAutoSave(300e9);
  fMitGPTree.tree_->SetDirectory(fOutputFile);
  AddOutput(fMitGPTree.tree_);

  // to the output module
  fMitGPTree.CreateTree(0);
  fMitGPTree.tree_->SetAutoSave(300e9);
  AddOutput(fMitGPTree.tree_);
}

//--------------------------------------------------------------------------------------------------
void BoostedVTreeWriter::SlaveTerminate()
{
  // Say how many events were analyzed
  printf("\n BoostedVTreeWriter::Terminate - Events analyzed: %d\n\n",fNAnalyzed);

  // Save the ntuple file
  fOutputFile->WriteTObject(fMitGPTree.tree_,fMitGPTree.tree_->GetName());
}

//--------------------------------------------------------------------------------------------------
float BoostedVTreeWriter::GetTau(fastjet::PseudoJet &iJet,int iN, float iKappa)
{
  // Calculate the tau variable for the given pseudojet (reproducing 1st paper)
  fastjet::contrib::Nsubjettiness nSubNKT(iN,fastjet::contrib::Njettiness::kt_axes,
					  iKappa,fConeSize,std::numeric_limits<double>::max());
  //fastjet::contrib::Nsubjettiness nSubNKT(iN,fastjet::contrib::Njettiness::onepass_kt_axes,
  //                                        iKappa,fConeSize,fConeSize);
  
  return nSubNKT(iJet);
}

//--------------------------------------------------------------------------------------------------
void BoostedVTreeWriter::ProcessMc()
{
  // We only get here for MC, so no more checking - we will perfrom the analysis on generator level
  // and take care of filling the tree

  LoadEventObject(fMcPartsName,fMcParts);

  // Fill Fastjet with the MC particles (no pileup included here)

  int nParts = 0;
  int nLeptons = 0;
  std::vector<fastjet::PseudoJet> lFjParts;
  for (UInt_t i=0; i<fMcParts->GetEntries(); ++i) {
    const MCParticle *p = fMcParts->At(i);
    // Grab out the first two leptons (el, mu, tau) that come from the W
    if (p->Status() == 3 &&
    	(p->Is(MCParticle::kEl) || p->Is(MCParticle::kMu) || p->Is(MCParticle::kTau)) &&
	p->Mother()->Is(MCParticle::kW)) {
      if      (nLeptons == 0) {
    	fMitGPTree.genLep1_.SetPxPyPzE(p->Px(),p->Py(),p->Pz(),p->E());
    	fMitGPTree.genLep1Pid_ = p->PdgId();
      }
      else if (nLeptons == 1) {
    	fMitGPTree.genLep2_.SetPxPyPzE(p->Px(),p->Py(),p->Pz(),p->E());
    	fMitGPTree.genLep2Pid_ = p->PdgId();
      }
      nLeptons++;
    }
    // Fill all particles into fastjet
    if (p->Status() == 1) { // just the particles that go to the detector
      nParts++;
      // Push all particles that make detector entries into fastjet particle collection
      lFjParts.push_back(fastjet::PseudoJet(p->Px(),p->Py(),p->Pz(),p->E()));
      lFjParts.back().set_user_index(i);
    }
  }

  // ---- Fastjet is ready ----

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea *lClustering =
    new fastjet::ClusterSequenceArea(lFjParts,*fCAJetDef,*fAreaDefinition);

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(lClustering->inclusive_jets(0.0));
 
  // Attach jets 1 and 2
  fastjet::PseudoJet lJet1, lJet2;
  if (fPrune) {
    if (lOutJets.size() > 0)
      lJet1 = (*fPruner)(lOutJets[0]);
    if (lOutJets.size() > 1)
      lJet2 = (*fPruner)(lOutJets[1]);
  }
  else {
    if (lOutJets.size() > 0)
      lJet1 = lOutJets[0];
    if (lOutJets.size() > 1)
      lJet2 = lOutJets[1];
  }

  // ---- Fastjet is already done ----

  // Things to store before cuts
  fMitGPTree.nGenParts_ = nParts;
  fMitGPTree.numGenJets_ = lOutJets.size();

  // Skip event if no jets are found
  if (lOutJets.size() == 0) {
    if (lClustering)
      delete lClustering;
    return;
  }

  // Basic cut on the first jet
  if (lJet1.pt() < 100) {
    if (lClustering)
      delete lClustering;
    return;
  }
  
  // Fill up the tree
  fMitGPTree.genJet1_.SetPxPyPzE(lJet1.px(),lJet1.py(),lJet1.pz(),lJet1.e());
  fMitGPTree.genJet1NParts_ = lJet1.constituents().size();
  fMitGPTree.genJet1Pt_ = lJet1.pt();
  fMitGPTree.genJet1Eta_ = lJet1.eta();
  fMitGPTree.genJet1Phi_ = lJet1.phi();
  fMitGPTree.genJet1M_ = lJet1.m();
  fMitGPTree.genJet1Tau1_ = GetTau(lJet1,1,1);
  fMitGPTree.genJet1Tau2_ = GetTau(lJet1,2,1);
  fMitGPTree.genJet1Tau3_ = GetTau(lJet1,3,1);
  fMitGPTree.genJet1MinTrigDr_ = MinTriggerDeltaR(fMitGPTree.genJet1_);

  if (lOutJets.size() > 1) {
    fMitGPTree.genJet2_.SetPxPyPzE(lJet2.px(),lJet2.py(),lJet2.pz(),lJet2.e());
    fMitGPTree.genJet2NParts_ = lJet2.constituents().size();
    fMitGPTree.genJet2Pt_ = lJet2.pt();
    fMitGPTree.genJet2Eta_ = lJet2.eta();
    fMitGPTree.genJet2Phi_ = lJet2.phi();
    fMitGPTree.genJet2M_ = lJet2.m();
    fMitGPTree.genJet2Tau1_ = GetTau(lJet2,1,1);
    fMitGPTree.genJet2Tau2_ = GetTau(lJet2,2,1);
    fMitGPTree.genJet2Tau3_ = GetTau(lJet2,3,1);
    fMitGPTree.genJet2MinTrigDr_ = MinTriggerDeltaR(fMitGPTree.genJet2_);
  }

  if (lClustering)
    delete lClustering;

  return;
}

//--------------------------------------------------------------------------------------------------
void BoostedVTreeWriter::GetJetTriggerObjs()
{
  // Pick only the jet trigger objects out of all our trigger objects

  // Make sure to initialize
  fJetTriggerObjs.clear();

  // Prepare loop
  fTrigObjs = GetHLTObjects(fTriggerObjsName);
  if (! fTrigObjs)
    printf(" BoostedVTreeWriter::Process() - ERROR - TriggerObjectCol not found\n");
  else {
    // loop through the stored trigger objects and find corresponding trigger name
    for (UInt_t j=0; j<fTrigObjs->GetEntries();++j) {
      const TriggerObject *to = fTrigObjs->At(j);
      TString trName = to->TrigName();
      // default MonoJet
      if (trName.Contains("MonoCentralPFJet80_PFMETnoMu"))
	fJetTriggerObjs.push_back(to);
    }
  }

  return;
}

//--------------------------------------------------------------------------------------------------
Double_t BoostedVTreeWriter::MinTriggerDeltaR(LorentzVector jet)
{
  // Determine the minimum delta R between the given provided jet and the jet trigger objects. This
  // will allow us to perfrom a trigger matching analysis.

  Double_t dR = 999;
  for (UInt_t i=0; i<fJetTriggerObjs.size(); i++) {
    Double_t dRTest = MathUtils::DeltaR(jet,*fJetTriggerObjs[i]);
    if (dRTest<dR)
      dR = dRTest;
  }

  return dR;
}

//--------------------------------------------------------------------------------------------------
bool BoostedVTreeWriter::IsTightMuon(const Muon *muon)
{
  return(((muon->HasGlobalTrk()                                   &&
           muon->GlobalTrk()->Chi2()/muon->GlobalTrk()->Ndof() < 10 &&
           (muon->NSegments() > 1 || muon->NMatches() > 1)          &&
           muon->NValidHits() > 0                                   ) ||
          muon->IsTrackerMuon()                                          ) &&
         (muon->BestTrk() != 0 && muon->BestTrk()->NHits() > 10 &&
          (muon->NSegments() > 1 || muon->NMatches() > 1)       &&
          muon->BestTrk()->NPixelHits() > 0                     )             );
}
