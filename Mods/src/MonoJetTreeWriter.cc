#include <TSystem.h>
#include <TFile.h>

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/GenericParticle.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PFMetCorrectionTools.h"

#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"

using namespace mithep;

ClassImp(mithep::MonoJetTreeWriter)

//--------------------------------------------------------------------------------------------------
MonoJetTreeWriter::MonoJetTreeWriter(const char *name, const char *title) :
  BaseMod                 (name,title),
  fEvtSelDataName         (Names::gkEvtSelDataBrn),
  fRawMetName             ("PFMet"),
  fMetName                ("PFMet"),
  fPhotonsName            (Names::gkPhotonBrn),
  fElectronsName          (Names::gkElectronBrn),
  fMuonsName              (Names::gkMuonBrn),
  fTausName               (Names::gkPFTauBrn),
  fJetsName               (Names::gkPFJetBrn),
  fRawJetsName            (Names::gkPFJetBrn),
  fLeptonsName            (ModNames::gkMergedLeptonsName),
  fPFCandidatesName       (Names::gkPFCandidatesBrn),
  fVertexName             (ModNames::gkGoodVertexesName),
  fSuperClustersName      ("PFSuperClusters"),
  fTracksName             (Names::gkTrackBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),
  fBeamspotName           (Names::gkBeamSpotBrn),
  fMCEvInfoName           (Names::gkMCEvtInfoBrn),
  fMCPartName             (Names::gkMCPartBrn),
  fTriggerObjectsName     ("MyHltPhotObjs"),
  fPFNoPileUpName         ("pfnopileupcands"),
  fPFPileUpName           ("pfpileupcands"),
  fIsData                 (false),
  fMetFromBranch          (kTRUE),
  fPhotonsFromBranch      (kTRUE),
  fElectronsFromBranch    (kTRUE),
  fMuonsFromBranch        (kTRUE),
  fTausFromBranch         (kTRUE),
  fJetsFromBranch         (kTRUE),
  fPFCandidatesFromBranch (kTRUE),
  fPVFromBranch           (kTRUE),
  fQGTaggerCHS            (kTRUE),
  // -------------------------
  fRawMet                 (0),
  fMet                    (0),
  fPhotons                (0),
  fElectrons              (0),
  fMuons                  (0),
  fPFTaus                 (0),
  fJets                   (0),
  fRawJets                (0),
  fTrigObj                (0),
  fPFCandidates           (0),
  fTracks                 (0),
  fPV                     (0),
  fBeamspot               (0),
  fMCEventInfo            (0),
  fPileUp                 (0),
  fPileUpDen              (0),
  fSuperClusters          (0),
  fParticles              (0),
  fEvtSelData             (0),
  // -------------------------
  fDecay                  (0),
  fFillNtupleType         (0),
  fNEventsSelected        (0),
  fOutputFile             (0)

{
  // WARNING, defining the object here invalidates the call of the setter for the CHS flag
  qgTagger = new QGTagger(fQGTaggerCHS);

  // Constructor
}

MonoJetTreeWriter::~MonoJetTreeWriter()
{
  // Destructor
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void MonoJetTreeWriter::SlaveTerminate()
{
  fOutputFile->WriteTObject(fMitGPTree.tree_,fMitGPTree.tree_->GetName());
  std::cout << "Processed events on MonoJetTreeWriter: " << fNEventsSelected << std::endl;
  //  delete fMVAMet;
}


//--------------------------------------------------------------------------------------------------
void MonoJetTreeWriter::Process()
{
  // Process entries of the tree.
  LoadEventObject(fBeamspotName,      fBeamspot);
  LoadEventObject(fEvtSelDataName,    fEvtSelData,    true);
  LoadEventObject(fPileUpDenName,     fPileUpDen,     true);
  MetOArr *GenMet = 0;
  MCParticleOArr *GenLeptons = 0;
  if (!fIsData) {
    LoadBranch(fMCEvInfoName);
    LoadBranch(fPileUpName);
    LoadEventObject(fMCPartName,      fParticles);
    GenMet = GetObject<MetOArr>(ModNames::gkMCMETName);
    GenLeptons = GetObject<MCParticleOArr>(ModNames::gkMCLeptonsName);
  }

  LoadEventObject(fRawMetName,        fRawMet,        true);
  LoadEventObject(fMetName,           fMet,           fMetFromBranch);
  LoadEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch);
  LoadEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  LoadEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  LoadEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  LoadEventObject(fJetsName,          fJets,          fJetsFromBranch);
  LoadEventObject(fRawJetsName,       fRawJets,       false);
  LoadEventObject(fPVName,            fPV,            fPVFromBranch);

  LoadEventObject(fPFCandidatesName,  fPFCandidates,  fPFCandidatesFromBranch);
  LoadEventObject(fSuperClustersName, fSuperClusters);
  LoadEventObject(fTracksName,        fTracks,        true);

  ParticleOArr    *leptons  = GetObject<ParticleOArr>(ModNames::gkMergedLeptonsName);
  const VertexCol *vertices = GetObject<VertexOArr>(fVertexName);

  const PFCandidateCol *fPFNoPileUpCands = GetObject<PFCandidateCol>(fPFNoPileUpName);    
  const PFCandidateCol *fPFPileUpCands = GetObject<PFCandidateCol>(fPFPileUpName);

  fNEventsSelected++;

  // initialize the tree variables
  fMitGPTree.InitVariables();

  // MC INFORMATION

  Double_t Q = 0.0;
  Int_t    id1 = 0;
  Double_t x1 = 0.0;
  Double_t pdf1 = 0.0;
  Int_t    id2 = 0;
  Double_t x2 = 0.0;
  Double_t pdf2 = 0.0;
  Int_t    processId = 0;
  if (! fIsData) {
    Q         = fMCEventInfo->Scale();
    id1       = fMCEventInfo->Id1();
    x1        = fMCEventInfo->X1();
    pdf1      = fMCEventInfo->Pdf1();
    id2       = fMCEventInfo->Id2();
    x2        = fMCEventInfo->X2();
    pdf2      = fMCEventInfo->Pdf2();
    processId = fMCEventInfo->ProcessId();

    // the Z
    for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
      if (fParticles->At(i)->Status()==3 and fParticles->At(i)->Is(MCParticle::kZ)) 
	fMitGPTree.genZ_ = fParticles->At(i)->Mom();
      if (fParticles->At(i)->Status()==3 and fParticles->At(i)->Is(MCParticle::kH)) 
	fMitGPTree.genH_ = fParticles->At(i)->Mom();

    }

    // the muons
    if(GenLeptons){
      if(GenLeptons->GetEntries() >= 1) {
	if (GenLeptons->At(0)->AbsPdgId()==13) fMitGPTree.genMuon1_ = GenLeptons->At(0)->Mom();
      }
      if(GenLeptons->GetEntries() >= 2) {
	if (GenLeptons->At(1)->AbsPdgId()==13) fMitGPTree.genMuon2_ = GenLeptons->At(1)->Mom();
      }
    }
  
    // the gen met
    fMitGPTree.genmet_ = GenMet->At(0)->Pt();
    fMitGPTree.genmetPhi_ = GenMet->At(0)->Phi();
  }

  fMitGPTree.Q_ = Q;
  fMitGPTree.id1_ = id1;
  fMitGPTree.x1_ = x1;
  fMitGPTree.pdf1_ = pdf1;
  fMitGPTree.id2_ = id2;
  fMitGPTree.x2_ = x2;
  fMitGPTree.pdf2_ = pdf2;
  fMitGPTree.processId_ = processId;

  // EVTSELDATA

  fMitGPTree.metFiltersWord_ = fEvtSelData->metFiltersWord();

  // PILEUP RELATED

  if (! fIsData) {
    // loop over the pileup summary info and grab what you need
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing() ==  0)
        fMitGPTree.npu_ = puinfo->GetPU_NumMean();
      if (puinfo->GetBunchCrossing() ==  1)
        fMitGPTree.npuPlusOne_ = puinfo->GetPU_NumInteractions();
      if (puinfo->GetBunchCrossing() == -1)
        fMitGPTree.npuMinusOne_ = puinfo->GetPU_NumInteractions();
    }
  }

  // TRIGGER

  fMitGPTree.trigger_ = 0;
  fTrigObj = GetHLTObjects(fTriggerObjectsName);

  if (! fTrigObj)
    printf("MonoJetTreeWriter::TriggerObjectCol not found\n");

  else {
    // loop through the stored trigger objects and find corresponding trigger name
    for (UInt_t i=0;i<fTrigObj->GetEntries();++i) {
      const TriggerObject *to = fTrigObj->At(i);
      TString trName = to->TrigName();
      // default MonoJet
      if (trName.Contains("MonoCentralPFJet80_PFMETnoMu"))
        fMitGPTree.trigger_ |= 1 << 0;
      if (trName.Contains("HLT_MET120_HBHENoiseCleaned_v") )
      fMitGPTree.trigger_ |= 1 << 1;
      // default VBF
      if (trName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"))
        fMitGPTree.trigger_ |= 1 << 2;
      // parked VBF, B-D. FIXME
      if (trName.Contains("HLT_DiJet30_MJJ700_AllJets_DEta3p5_VBF_v1") ||
          trName.Contains("HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v5"))
        fMitGPTree.trigger_ |= 1 << 3;
      // default single muon
      if (trName.Contains("HLT_IsoMu24_eta2p1_v"))
        fMitGPTree.trigger_ |= 1 << 4;
    }
  }

  fMitGPTree.run_ = GetEventHeader()->RunNum();
  fMitGPTree.lumi_ = GetEventHeader()->LumiSec();
  fMitGPTree.event_ = GetEventHeader()->EvtNum();
  fMitGPTree.nvtx_ = fPV->GetEntries();
  fMitGPTree.scale1fb_ = 1000.0;
  fMitGPTree.cuts_ = MitGPTree::DiLepton;

  if (fDecay == 0)
    fMitGPTree.dstype_ = MitGPTree::data;
  else
    fMitGPTree.dstype_ = MitGPTree::other;

  // MET BASICS

  fMitGPTree.metRaw_        = fRawMet->At(0)->Pt();
  fMitGPTree.metRawPhi_     = fRawMet->At(0)->Phi();
  fMitGPTree.metSig_        = fRawMet->At(0)->PFMetSig(); // RawPFMet for sig. calculation
  fMitGPTree.met_           = fMet->At(0)->Pt();
  fMitGPTree.metPhi_        = fMet->At(0)->Phi();
  fMitGPTree.sumEt_         = fMet->At(0)->SumEt();

  // these will be recalculated removing one lepton (W - W boson)
  fMitGPTree.metRawCorW_    = fRawMet->At(0)->Pt();       // default value as RawMet
  fMitGPTree.metRawCorWPhi_ = fRawMet->At(0)->Phi();      // default value as RawMet
  fMitGPTree.metCorW_       = fMet->At(0)->Pt();          // default value as Met
  fMitGPTree.metCorWPhi_    = fMet->At(0)->Phi();         // default value as Met
  // these will be recalculated removing two leptons (Z - Z boson)
  fMitGPTree.metRawCorZ_    = fRawMet->At(0)->Pt();       // default value as RawMet
  fMitGPTree.metRawCorZPhi_ = fRawMet->At(0)->Phi();      // default value as RawMet
  fMitGPTree.metCorZ_       = fMet->At(0)->Pt();          // default value as Met
  fMitGPTree.metCorZPhi_    = fMet->At(0)->Phi();         // default value as Met

  // LEPTONS

  fMitGPTree.nlep_ = leptons->GetEntries();
  if (leptons->GetEntries() >= 1) {         // loop over all leptons
    const Particle *lep = leptons->At(0);
    
    fMitGPTree.lep1_ = lep->Mom();
    if      (lep->ObjType() == kMuon) {
      fMitGPTree.lid1_ = 13;
      const Muon* mu = dynamic_cast<const Muon*>(lep);
      fMitGPTree.lep1IsTightMuon_ = IsTightMuon(mu);
      fMitGPTree.lep1PtErr_ = mu->BestTrk()->PtErr()/mu->BestTrk()->Pt();
      double totalIso =  IsolationTools::BetaMwithPUCorrection(fPFNoPileUpCands, fPFPileUpCands, mu, 0.4);
      fMitGPTree.lep1IsIsolated_ = totalIso < (mu->Pt()*0.2);
    }
    else if (lep->ObjType() == kElectron)
      fMitGPTree.lid1_ = 11;
    else
      assert(0);  // cannot happen: leptons are only muons and electrons

    if (lep->Charge() < 0)
      fMitGPTree.lid1_ = -1 * fMitGPTree.lid1_;

    // If the event contains at least 1 lepton correct the MET using the highest pt one
    CorrectMet(fMitGPTree.met_,    fMitGPTree.metPhi_,leptons->At(0),0,
	       fMitGPTree.metCorW_,fMitGPTree.metCorWPhi_);
    CorrectMet(fMitGPTree.metRaw_, fMitGPTree.metRawPhi_,leptons->At(0),0,
	       fMitGPTree.metRawCorW_,fMitGPTree.metRawCorWPhi_);
  }

  if (leptons->GetEntries() >= 2) {
    const Particle *lep = leptons->At(1);

    fMitGPTree.lep2_ = lep->Mom();
    if     (lep->ObjType() == kMuon) {
      fMitGPTree.lid2_ = 13;
      const Muon* mu = dynamic_cast<const Muon*>(lep);
      fMitGPTree.lep2IsTightMuon_ = IsTightMuon(mu);
      fMitGPTree.lep2PtErr_ = mu->BestTrk()->PtErr()/mu->BestTrk()->Pt();
      double totalIso =  IsolationTools::BetaMwithPUCorrection(fPFNoPileUpCands, fPFPileUpCands, mu, 0.4);
      fMitGPTree.lep2IsIsolated_ = totalIso < (mu->Pt()*0.2);
    }
    else if (lep->ObjType() == kElectron)
      fMitGPTree.lid2_ = 11;
    else
      assert(0);

    if (lep->Charge() < 0)
      fMitGPTree.lid2_ = -1 * fMitGPTree.lid2_;

    // If the event contains at least 2 leptons correct the MET using the 2 highest pt ones
    CorrectMet(fMitGPTree.met_,   fMitGPTree.metPhi_,leptons->At(0),leptons->At(1),
	       fMitGPTree.metCorZ_,fMitGPTree.metCorZPhi_);
    CorrectMet(fMitGPTree.metRaw_,fMitGPTree.metRawPhi_,leptons->At(0),leptons->At(1),
	       fMitGPTree.metRawCorZ_,fMitGPTree.metRawCorZPhi_);
  }

  if (leptons->GetEntries() >= 3) {
    const Particle *lep = leptons->At(2);

    fMitGPTree.lep3_ = lep->Mom();
    if      (lep->ObjType() == kMuon) {
      fMitGPTree.lid3_ = 13;
      fMitGPTree.lep3IsTightMuon_ = IsTightMuon(dynamic_cast<const Muon*>(leptons->At(2)));
    }
    else if (lep->ObjType() == kElectron)
      fMitGPTree.lid3_ = 11;
    else
      assert(0);

    if (lep->Charge() < 0)
      fMitGPTree.lid3_ = -1 * fMitGPTree.lid3_;

    // If the event contains more than 2 we have no further assumption ( or should we? WZ ;-) )

  }

  // PHOTON(S)

  fMitGPTree.nphotons_ = fPhotons->GetEntries();
  if (fPhotons->GetEntries() >= 1) {
    const Photon *photon = fPhotons->At(0);
    fMitGPTree.pho1_ = photon->Mom();
    if (fPhotons->GetEntries() >= 2) {
      photon = fPhotons->At(1);
      fMitGPTree.pho2_ = photon->Mom();
    }
  }

  // TAUS

  fMitGPTree.ntaus_ = fPFTaus->GetEntries();
  if (fPFTaus->GetEntries() >= 1) {
    const PFTau *tau = fPFTaus->At(0);
    fMitGPTree.tau1_ = tau->Mom();
    if (fPFTaus->GetEntries() >= 2) {
      tau = fPFTaus->At(1);
      fMitGPTree.tau2_ = tau->Mom();
    }
  }

  // JETS

  const TriggerObject *trigObj1 = 0;
  const TriggerObject *trigObj2 = 0;

  // for mvamet (?)
  PFJetOArr *pfJets = new PFJetOArr;
  pfJets->SetOwner(kTRUE);
  // pfJets->SetName(fPFJetsName);
  for (UInt_t i=0; i<fRawJets->GetEntries(); ++i) {
    const Jet *inJet = fRawJets->At(i);

    // copy input jet, using special function to a deep copy (and own it)
    Jet *jet = inJet->MakeCopy();
    pfJets->AddOwned(dynamic_cast<PFJet*>(jet));

  }

  fMitGPTree.njets_ = fJets->GetEntries();
  fMitGPTree.noiseCleaning_ = 0;
  qgTagger->SetRhoIso(fPileUpDen->At(0)->RhoRandomLowEta());

  if (fJets->GetEntries() >= 1) {
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(0));

    fMitGPTree.jet1_ = jet->Mom();

    fMitGPTree.jet1CHF_ = jet->ChargedHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet1NHF_ = jet->NeutralHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet1NEMF_ = jet->NeutralEmEnergy()/jet->RawMom().E();
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet1CHF_>0.2) << 0;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet1NHF_>0.7) << 1;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet1NEMF_>0.7) << 2;
    fMitGPTree.jet1Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();

    qgTagger->CalculateVariables(jet, vertices);
    fMitGPTree.jet1QGtag_ = qgTagger->QGValue();

    // variables used for the QG retraining
    fMitGPTree.jet1QGRho_   = fPileUpDen->At(0)->RhoRandomLowEta();
    fMitGPTree.jet1QGPtD_   = qgTagger->GetPtD();
    fMitGPTree.jet1QGAxis1_ = qgTagger->GetAxis1();
    fMitGPTree.jet1QGAxis2_ = qgTagger->GetAxis2();
    fMitGPTree.jet1QGMult_  = qgTagger->GetMult();

    // matching in MC
    if (! fIsData) {
      double minPartonicDR = 0.8;
      UInt_t partonId = 0;
      for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
        const MCParticle *p = fParticles->At(i);
        if (p->Status()!=3)
          continue;
        if (p->AbsPdgId()>5 and p->AbsPdgId()!=21)
          continue; //1, 2, 3, 4, 5, 21
        if (MathUtils::DeltaR(*p,*jet)< minPartonicDR) {
          minPartonicDR = MathUtils::DeltaR(*p,*jet);
          partonId = p->AbsPdgId();
        }
      }
      fMitGPTree.jet1PartonId_ = partonId;
    }

    // trigger matching
    for (UInt_t i=0; i<fTrigObj->GetEntries(); ++i){
      const TriggerObject *trigobj=fTrigObj->At(i);
      TString trigName = trigobj->TrigName();
      if (trigobj->IsHLT() ){
        if (trigName.Contains("MonoCentralPFJet80_PFMETnoMu")) {
          bool match = true;
          if (trigobj->Pt() < 80)
            match = false;
          if (trigobj->Type() != 85)
            match = false;
          if (MathUtils::DeltaR(trigobj->Phi(),trigobj->Eta(),jet->Phi(),jet->Eta()) > 0.5)
            match = false;
          if (match)
            fMitGPTree.HLTmatch_ |= 1 << 0;
        }
        if (trigName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v")) {
          bool match = true;
          if (trigobj->Pt() < 40)
            match = false;
          if (trigobj->Type() != 85)
            match = false;
          if (MathUtils::DeltaR(trigobj->Phi(), trigobj->Eta(), jet->Phi(), jet->Eta()) > 0.5)
            match = false;
          if (match)
	    trigObj1 = trigobj;
        }
      }
    }
  }

  if (fJets->GetEntries() >= 2) {
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(1));
    fMitGPTree.jet2_ = jet->Mom();
    fMitGPTree.jet2CHF_  = jet->ChargedHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet2NHF_  = jet->NeutralHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet2NEMF_  = jet->NeutralEmEnergy()/jet->RawMom().E();
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet2CHF_>0.2) << 3;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet2NHF_>0.7) << 4;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet2NEMF_>0.7) << 5;
    fMitGPTree.jet2Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
    qgTagger->CalculateVariables(jet, vertices);
    fMitGPTree.jet2QGtag_ = qgTagger->QGValue();

    // trigger matching
    for (UInt_t i=0; i<fTrigObj->GetEntries(); ++i) {
      const TriggerObject *trigobj = fTrigObj->At(i);
      TString trigName = trigobj->TrigName();
      if (trigobj->IsHLT() ) {
        if (trigName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v")) {
          bool match = true;
          if (trigobj->Pt() < 40)
            match = false;
          if (trigobj->Type() != 85)
            match = false;
          if (MathUtils::DeltaR(trigobj->Phi(),trigobj->Eta(),jet->Phi(),jet->Eta()) > 0.5)
            match = false;
          if (match)
	    trigObj2 = trigobj;
        }
      }
    }
  }
  if (trigObj1 && trigObj2) {
    double diJetMass = (trigObj1->Mom()+trigObj2->Mom()).M();
    double deltaEta  = fabs(trigObj1->Eta()-trigObj2->Eta());
    if ((diJetMass > 800) && (deltaEta > 3.5) && (trigObj1->Eta()*trigObj2->Eta() < 0)) {
      fMitGPTree.HLTmatch_ |= 1 << 2;
    }
  }
  if (fJets->GetEntries() >= 3) {
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(2));
    fMitGPTree.jet3_     = jet->Mom();
    fMitGPTree.jet3CHF_  = jet->ChargedHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet3NHF_  = jet->NeutralHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet3NEMF_ = jet->NeutralEmEnergy()/jet->RawMom().E();
    fMitGPTree.jet3Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }
  if (fJets->GetEntries() >= 4) {
    const Jet *jet = fJets->At(3);
    fMitGPTree.jet4_     = jet->Mom();
    fMitGPTree.jet4Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }

  // TRACKS

  fMitGPTree.ntracks_ = 0;
  for (UInt_t i=0; i<fTracks->GetEntries(); i++) {
    const mithep::Track* pTrack = fTracks->At(i);
    if (pTrack->Pt() <= 15)
      continue;

    Bool_t isLepton = kFALSE;
    for (UInt_t j=0; j<leptons->GetEntries(); j++) {
      const Particle *lep = leptons->At(j);
      if (MathUtils::DeltaR(pTrack->Mom(), lep->Mom()) < 0.05) {
        isLepton = kTRUE;
        break;
      }
    }
    if (isLepton == kTRUE)
      continue;

    GenericParticle *p = new GenericParticle(pTrack->Px(), pTrack->Py(), pTrack->Pz(),
                                             pTrack->P(), pTrack->Charge());
    if (fMitGPTree.ntracks_ == 0)
      fMitGPTree.track1_ = p->Mom();
    if (fMitGPTree.ntracks_ == 1)
      fMitGPTree.track2_ = p->Mom();
    if (fMitGPTree.ntracks_ == 2)
      fMitGPTree.track3_ = p->Mom();
    delete p;
    fMitGPTree.ntracks_++;
  }

  // Finally fill the tree
  fMitGPTree.tree_->Fill();

  delete pfJets;

  return;
}

//--------------------------------------------------------------------------------------------------
void MonoJetTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  if (! fIsData) {
    ReqEventObject(fMCPartName,      fParticles,     true);
    ReqBranch(fPileUpName,           fPileUp);
    ReqBranch(fMCEvInfoName,         fMCEventInfo);
  }
  ReqEventObject(fBeamspotName,      fBeamspot,      true);
  ReqEventObject(fPileUpDenName,     fPileUpDen,     true);
  ReqEventObject(fPVName,            fPV,            fPVFromBranch);
  ReqEventObject("EvtSelData",       fEvtSelData,    true);

  ReqEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch);
  ReqEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  ReqEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  ReqEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  ReqEventObject(fPFCandidatesName,  fPFCandidates,  fPFCandidatesFromBranch);
  ReqEventObject(fJetsName,          fJets,          fJetsFromBranch);
  ReqEventObject(fRawJetsName,       fRawJets,       false);
  ReqEventObject(fRawMetName,        fRawMet,        true);
  ReqEventObject(fMetName,           fMet,           fMetFromBranch);

  ReqEventObject(fSuperClustersName, fSuperClusters, true);
  ReqEventObject(fTracksName,        fTracks,        true);


  // This should be done in the run macro
  if (fIsData) {
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()));

    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data()));
  }
  else {
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()));
  }


  // Create Ntuple Tree
  fOutputFile = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  fMitGPTree.CreateTree(fFillNtupleType);
  fMitGPTree.tree_->SetAutoSave(300e9);
  fMitGPTree.tree_->SetDirectory(fOutputFile);
  AddOutput(fMitGPTree.tree_);
}

bool MonoJetTreeWriter::IsTightMuon(const Muon *muon)
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

void MonoJetTreeWriter::CorrectMet(const float met, const float metPhi,
				   const Particle *l1, const Particle *l2,
                                   float &newMet, float &newMetPhi)
{
  // inputs:  met, metPhi, l1,  [ l2  only used if pointer is non-zero ]
  // outputs: newMet, newMetPhi

  float newMetX;
  float newMetY;
  if (l2) {   // these are doubles on the right: need full calculation in one line
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px() + l2->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py() + l2->Mom().Py();
  }
  else {
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py();
  }
  
  newMet    = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);
}
