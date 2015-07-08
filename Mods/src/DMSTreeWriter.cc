#include <TSystem.h>
#include <TFile.h>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/GenericParticle.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/XlJet.h"
#include "MitAna/DataTree/interface/XlFatJet.h"
#include "MitAna/DataTree/interface/XsIsoParticle.h"
#include "MitMonoJet/Utils/interface/DiJetMVA.h"
#include "MitAna/DataTree/interface/XlSubJet.h"

#include "MitMonoJet/Mods/interface/DMSTreeWriter.h"

using namespace mithep;

ClassImp(mithep::DMSTreeWriter)

//--------------------------------------------------------------------------------------------------
DMSTreeWriter::DMSTreeWriter(const char *name, const char *title) :
  BaseMod                 (name,title),
  fEvtSelDataName         ("XlEvtSelData"),
  fRawMetName             ("PFMet"),
  fPhotonsName            ("XsPhotons"),
  fElectronsName          ("XsElectrons"),
  fMuonsName              ("XsMuons"),
  fTausName               ("XsTaus"),
  fJetsName               ("XlJets"),
  fFatJetsName            ("XlFatJets"),
  fSubJetsName            ("XlSubJets"),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),
  fMCEventInfoName        (Names::gkMCEvtInfoBrn),
  fMCParticlesName        (Names::gkMCPartBrn),
  fTriggerObjectsName     ("MyHltPhotObjs"),
  fIsData                 (kFALSE),
  fPhotonsFromBranch      (kTRUE),
  fElectronsFromBranch    (kTRUE),
  fMuonsFromBranch        (kTRUE),
  fTausFromBranch         (kTRUE),
  fJetsFromBranch         (kTRUE),
  fFatJetsFromBranch      (kTRUE),
  fSubJetsFromBranch      (kTRUE),
  fPVFromBranch           (kTRUE),
  // -------------------------
  fRawMet                 (0),
  fPhotons                (0),
  fElectrons              (0),
  fMuons                  (0),
  fPFTaus                 (0),
  fJets                   (0),
  fFatJets                (0),
  fSubJets                (0),
  fPV                     (0),
  fPileUp                 (0),
  fPileUpDen              (0),
  fMCEventInfo            (0),
  fMCParticles            (0),
  fEvtSelData             (0),
  fTrigObj                (0),
  fPUInputFileName        ("MyInputPUFile"),       
  fPUTargetFileName       ("MyTargetPUFile"),
  fPUInput                (0),       
  fPUTarget               (0),
  fPUWeight               (0),
  fDiJetMVA               (0),
  fJetUncertainties       (0),
  fFatJetUncertainties    (0),
   // -------------------------
  fOutputFile             (0)

{
  // Constructor
}

DMSTreeWriter::~DMSTreeWriter()
{
  // Destructor
  delete fDiJetMVA;
  delete fJetUncertainties;
  delete fFatJetUncertainties;
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::SlaveTerminate()
{
  fOutputFile->WriteTObject(fMitDMSTree.tree_,fMitDMSTree.tree_->GetName());
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::Process()
{
  // Process entries of the tree.
  LoadEventObject(fEvtSelDataName,    fEvtSelData,    true);
  LoadEventObject(fPileUpDenName,     fPileUpDen,     true);
  if (!fIsData) {
    LoadBranch(fPileUpName);
    LoadBranch(fMCEventInfoName);
    LoadBranch(fMCParticlesName);
  }
  LoadEventObject(fTriggerObjectsName,fTrigObj,       true);

  LoadEventObject(fRawMetName,        fRawMet,        true);
  LoadEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch);
  LoadEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  LoadEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  LoadEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  LoadEventObject(fJetsName,          fJets,          fJetsFromBranch);
  LoadEventObject(fFatJetsName,       fFatJets,       fJetsFromBranch);
  LoadEventObject(fSubJetsName,       fSubJets,       fJetsFromBranch);
  LoadEventObject(fPVName,            fPV,            fPVFromBranch);

  // initialize the tree variables
  fMitDMSTree.InitVariables();
  
  // EVTSELDATA  
  fMitDMSTree.metFiltersWord_ = fEvtSelData->metFiltersWord();
  fMitDMSTree.trigger_        = fEvtSelData->HLTWord();
  fMitDMSTree.preselWord_     = fEvtSelData->preselWord();

  // PILEUP RELATED

  if (! fIsData) {
    // loop over the pileup summary info and grab what you need
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing() ==  0)
        fMitDMSTree.npu_ = puinfo->GetPU_NumMean();
      if (puinfo->GetBunchCrossing() ==  1)
        fMitDMSTree.npuPlusOne_ = puinfo->GetPU_NumInteractions();
      if (puinfo->GetBunchCrossing() == -1)
        fMitDMSTree.npuMinusOne_ = puinfo->GetPU_NumInteractions();
    }
  }

  // EVENT

  fMitDMSTree.run_   = GetEventHeader()->RunNum();
  fMitDMSTree.lumi_  = GetEventHeader()->LumiSec();
  fMitDMSTree.event_ = GetEventHeader()->EvtNum();
  fMitDMSTree.nvtx_  = fPV->GetEntries();

  // HLT CHECK
  if (! fTrigObj)
    printf("MonoJetTreeWriter::TriggerObjectCol not found\n");

  // MET BASICS

  fMitDMSTree.metRaw_        = fRawMet->At(0)->Pt();
  fMitDMSTree.metRawPhi_     = fRawMet->At(0)->Phi();
  fMitDMSTree.met_           = fMitDMSTree.metRaw_;
  fMitDMSTree.metPhi_        = fMitDMSTree.metRawPhi_;
      
  // LEPTONS (MU and ELE), save tight id for further studies
  // Also perform met/mt/mll computation according to the relevant case (muons only)
  fMitDMSTree.nele_ = fElectrons->GetEntries();
  fMitDMSTree.nlep_ = fMuons->GetEntries();
  if (fMuons->GetEntries() > 1) {
    //mu mu
    fMitDMSTree.lep1_ = fMuons->At(0)->Mom();
    fMitDMSTree.lid1_ = fMuons->At(0)->Charge()*13;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eTightMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*10;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eIsoMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*100;

    fMitDMSTree.lep2_ = fMuons->At(1)->Mom();
    fMitDMSTree.lid2_ = fMuons->At(1)->Charge()*13;
    if (fMuons->At(1)->ParticleId() == XsIsoParticle::EParticleId::eTightMuon)
      fMitDMSTree.lid2_ = fMitDMSTree.lid2_*10;
    if (fMuons->At(1)->ParticleId() == XsIsoParticle::EParticleId::eIsoMuon)
      fMitDMSTree.lid2_ = fMitDMSTree.lid2_*100;

    CorrectMet(fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_,fMitDMSTree.lep1_,fMitDMSTree.lep2_,
               fMitDMSTree.met_,fMitDMSTree.metPhi_);    
    fMitDMSTree.mll_ = (fMitDMSTree.lep1_ + fMitDMSTree.lep2_).M();
  }
  else if (fMuons->GetEntries() == 1) {
    fMitDMSTree.lep1_ = fMuons->At(0)->Mom();
    fMitDMSTree.lid1_ = fMuons->At(0)->Charge()*13;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eTightMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*10;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eIsoMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*100;

    CorrectMet(fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_,fMitDMSTree.lep1_,
               fMitDMSTree.met_,fMitDMSTree.metPhi_);    
    fMitDMSTree.mt_ = GetMt(fMitDMSTree.lep1_,fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_);
  }
  // Make Muon-HLT matching 
  if (fMitDMSTree.nlep_ > 0 && fMitDMSTree.lid1_ >= 13
   && IsHLTMatched(fMitDMSTree.lep1_, 
                   TriggerObject::TriggerMuon,
                   0.3, 
                   24., 2.1))
    fMitDMSTree.HLTmatch_ |= MitDMSTree::MuonMatch;         

  // TAUS

  fMitDMSTree.ntaus_ = fPFTaus->GetEntries();
  if (fPFTaus->GetEntries() >= 1) {
    fMitDMSTree.tau1_ = fPFTaus->At(0)->Mom();
  }

  // PHOTON(S)

  fMitDMSTree.nphotons_ = fPhotons->GetEntries();
  if (fPhotons->GetEntries() >= 1) {
    const XsIsoParticle *photon = fPhotons->At(0);
    fMitDMSTree.pho1_ = photon->Mom();
    // Make Photon-HLT matching 
    if (IsHLTMatched(fMitDMSTree.pho1_, 
                     TriggerObject::TriggerPhoton,
                     0.3, 
                     130.))
      fMitDMSTree.HLTmatch_ |= MitDMSTree::PhotonMatch;         
    // Correct MET (only if no leptons and high pt photon!) and correct for Fprint
    if (fMitDMSTree.nele_ == 0 && fMitDMSTree.nlep_ == 0 && photon->Pt() > 150.) {
      CorrectMet(fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_,fMitDMSTree.pho1_,
                 fMitDMSTree.met_,fMitDMSTree.metPhi_);    
      CorrectMet(fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_,fMitDMSTree.pho1_,
                 fMitDMSTree.metFprint_,fMitDMSTree.metFprintPhi_,true);
    }    
  }
  
  // JETS : careful since the hardest could overlap with the fat jets
  fMitDMSTree.rmvaval_ = -10.;
  float resolvedVars[7];
  UInt_t resolvedIndexOne = 0;
  UInt_t resolvedIndexTwo = 0;
  fMitDMSTree.njets_ = 0;
  fMitDMSTree.njetsUp_ = 0;
  fMitDMSTree.njetsDown_ = 0;
  fMitDMSTree.nbjets_ = 0;
  for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
    
    const XlJet *jet = fJets->At(i);
    // Perform jet cleaning according to lepton counting
    if (fMitDMSTree.nlep_ > 0 
     && MathUtils::DeltaR(fMitDMSTree.lep1_,jet->Mom()) < 0.5)
      continue;
    if (fMitDMSTree.nlep_ > 1 
     && MathUtils::DeltaR(fMitDMSTree.lep2_,jet->Mom()) < 0.5)
      continue;

    // Get the jet uncertainty
    fJetUncertainties->setJetPt(jet->Pt());
    fJetUncertainties->setJetEta(jet->Eta());
    float jetUnc = fJetUncertainties->getUncertainty(true);
    
    if (fMitDMSTree.njets_ == 0) {
      fMitDMSTree.jet1_        = jet->Mom();
      fMitDMSTree.jet1Unc_     = jetUnc;
      fMitDMSTree.jet1CHF_     = jet->ChargedHadronEnergy()/jet->RawMom().E();
      fMitDMSTree.jet1NHF_     = jet->NeutralHadronEnergy()/jet->RawMom().E();
      fMitDMSTree.jet1NEMF_    = jet->NeutralEmEnergy()/jet->RawMom().E();
      // Make Jet-HLT matching: this is used for preselection
      if (IsHLTMatched(fMitDMSTree.jet1_, 
                       TriggerObject::TriggerJet,
                       0.5, 
                       80., 2.4))
        fMitDMSTree.HLTmatch_ |= MitDMSTree::JetMatch;         

      // Angular info
      fMitDMSTree.jet1metDphi_ = MathUtils::DeltaPhi(jet->Phi(),(double)fMitDMSTree.metPhi_);
      fMitDMSTree.jet1jet2Dphi_ = GetJetJetsDphi(fMitDMSTree.jet1_);
    }
    if (fMitDMSTree.njets_ == 1)
      fMitDMSTree.jet2_        = jet->Mom();
    if (fMitDMSTree.njets_ == 2)
      fMitDMSTree.jet3_        = jet->Mom();
    if (fMitDMSTree.njets_ == 3)
      fMitDMSTree.jet4_        = jet->Mom();
    if (fMitDMSTree.njets_ == 4)
      fMitDMSTree.jet5_        = jet->Mom();

    // Increment the cleaned jet counter and the relative syst counters
    fMitDMSTree.njets_++;
    if (jet->Pt()*(1.+jetUnc) > 30.) fMitDMSTree.njetsUp_++;
    if (jet->Pt()*(1.-jetUnc) > 30.) fMitDMSTree.njetsDown_++;

    // Start resolved section
    if (fJets->GetEntries() > i && (fMitDMSTree.preselWord_ & MitDMSTree::Resolved)) {
      if (jet->CombinedSecondaryVertexBJetTagsDisc() > 0.679)
        continue;
      for (UInt_t j = i+1; j < fJets->GetEntries(); ++j) {
        const XlJet *jetTwo = fJets->At(j);
        if (jetTwo->CombinedSecondaryVertexBJetTagsDisc() > 0.679)
          continue;
        // Perform jet cleaning according to lepton counting
        if (fMitDMSTree.nlep_ > 0 
         && MathUtils::DeltaR(fMitDMSTree.lep1_,jetTwo->Mom()) < 0.5)
          continue;
        if (fMitDMSTree.nlep_ > 1 
         && MathUtils::DeltaR(fMitDMSTree.lep2_,jetTwo->Mom()) < 0.5)
          continue;
        // Check mass
        if ((jet->Mom() + jetTwo->Mom()).M() < 60.
          ||(jet->Mom() + jetTwo->Mom()).M() > 110.) 
          continue;
        Double_t thisMVAval = fDiJetMVA->MVAValue(
                              jet,jetTwo,
                              fPileUpDen->At(0)->RhoRandomLowEta(),
                              fMCParticles,kFALSE,resolvedVars);
        if (thisMVAval > fMitDMSTree.rmvaval_) {
          fMitDMSTree.rmvaval_ = thisMVAval;
          resolvedIndexOne = i;
          resolvedIndexTwo = j;
          fMitDMSTree.rptOverM_     = resolvedVars[0];       
          fMitDMSTree.rjet1_pullang_= resolvedVars[1];
          fMitDMSTree.rjet2_pullang_= resolvedVars[2];
          fMitDMSTree.rjet1_qgl_    = resolvedVars[3];
          fMitDMSTree.rjet2_qgl_    = resolvedVars[4];
          fMitDMSTree.rmdrop_       = resolvedVars[5];
        }
      } // endl loop on second jet
    } // end resolved mini-ana

    // Check that the jet is b-tagged (medium WP)
    float btag = jet->CombinedSecondaryVertexBJetTagsDisc();  
    if (btag < 0.679)
      continue;

    // Fill the information for the two hardest b-jets
    if (fMitDMSTree.nbjets_ == 0) {
      fMitDMSTree.bjet1_        = jet->Mom();
      fMitDMSTree.bjet1Btag_    = btag;
    }
    if (fMitDMSTree.nbjets_ == 1) {
      fMitDMSTree.bjet2_        = jet->Mom();
      fMitDMSTree.bjet2Btag_    = btag;
    }  

    // Increment the b-jet counter
    fMitDMSTree.nbjets_ ++;

  }

  // Fill revant information in case a good diJet pair is found
  if (fMitDMSTree.rmvaval_ >= -1.) {
    fMitDMSTree.rjet1_ = fJets->At(resolvedIndexOne)->Mom();
    fMitDMSTree.rjet2_ = fJets->At(resolvedIndexTwo)->Mom();
  }

  // FAT JETS  
  fMitDMSTree.nfjets_ = 0;
  for (UInt_t i = 0; i < fFatJets->GetEntries(); ++i) {

    const XlFatJet *fjet = fFatJets->At(i);    
    // Perform jet cleaning according to lepton counting
    if (fMitDMSTree.nlep_ > 0 
     && MathUtils::DeltaR(fMitDMSTree.lep1_,fFatJets->At(i)->Mom()) < 0.5)
      continue;
    if (fMitDMSTree.nlep_ > 1 
     && MathUtils::DeltaR(fMitDMSTree.lep2_,fFatJets->At(i)->Mom()) < 0.5)
      continue;

    // Get the jet uncertainty
    fFatJetUncertainties->setJetPt(fjet->Pt());
    fFatJetUncertainties->setJetEta(fjet->Eta());
    float jetUnc = fFatJetUncertainties->getUncertainty(true);
          
    if (fMitDMSTree.nfjets_ == 0) {
      fMitDMSTree.fjet1_       = fjet->Mom();
      fMitDMSTree.fjet1Unc_     = jetUnc;
      fMitDMSTree.fjet1CHF_     = fjet->ChargedHadronEnergy()/fjet->RawMom().E();
      fMitDMSTree.fjet1NHF_     = fjet->NeutralHadronEnergy()/fjet->RawMom().E();
      fMitDMSTree.fjet1NEMF_    = fjet->NeutralEmEnergy()/fjet->RawMom().E();
      fMitDMSTree.fjet1Btag_    = GetFatJetBtag(fMitDMSTree.fjet1_, 0.5);
      fMitDMSTree.fjet1Charge_  = fjet->Charge();
      fMitDMSTree.fjet1QGtag_   = fjet->QGTag();
      fMitDMSTree.fjet1Tau1_    = fjet->Tau1();
      fMitDMSTree.fjet1Tau2_    = fjet->Tau2();
      fMitDMSTree.fjet1Tau3_    = fjet->Tau3();
      fMitDMSTree.fjet1C2b0_    = fjet->C2b0();
      fMitDMSTree.fjet1C2b0p2_      = fjet->C2b0p2();      
      fMitDMSTree.fjet1C2b0p5_      = fjet->C2b0p5();      
      fMitDMSTree.fjet1C2b1_        = fjet->C2b1();        
      fMitDMSTree.fjet1C2b2_        = fjet->C2b2();        
      fMitDMSTree.fjet1QJetVol_     = fjet->QJetVol();     
      fMitDMSTree.fjet1MassSDbm1_   = fjet->MassSDbm1();   
      fMitDMSTree.fjet1MassSDb0_    = fjet->MassSDb0();    
      fMitDMSTree.fjet1MassSDb1_    = fjet->MassSDb1();    
      fMitDMSTree.fjet1MassSDb2_    = fjet->MassSDb2();    
      fMitDMSTree.fjet1MassPruned_  = fjet->MassPruned();  
      fMitDMSTree.fjet1MassFiltered_= fjet->MassFiltered();
      fMitDMSTree.fjet1MassTrimmed_ = fjet->MassTrimmed();
      fMitDMSTree.fjet1Pull_        = fjet->Pull();
      fMitDMSTree.fjet1PullAngle_   = fjet->PullAngle();
      if (!fIsData)  
        fMitDMSTree.fjet1PartonId_  = JetPartonMatch(fMitDMSTree.fjet1_, 0.7);  
    
      // Subjets info
      fMitDMSTree.fjet1nsj_ = fjet->NSubJets();
      if (fMitDMSTree.fjet1nsj_ > 0) {
        fMitDMSTree.fjet1sj1_         = fjet->SubJet(0)->Mom();
        fMitDMSTree.fjet1QGtagSub1_   = fjet->SubJet(0)->QGTag();
        fMitDMSTree.fjet1QGPtDSub1_   = fjet->SubJet(0)->QGPtD();
        fMitDMSTree.fjet1QGAxis1Sub1_ = fjet->SubJet(0)->QGAxis1();
        fMitDMSTree.fjet1QGAxis2Sub1_ = fjet->SubJet(0)->QGAxis2();
        fMitDMSTree.fjet1QGMultSub1_  = fjet->SubJet(0)->QGMult();
      }
      if (fMitDMSTree.fjet1nsj_ > 1) {
        fMitDMSTree.fjet1sj2_         = fjet->SubJet(1)->Mom();
        fMitDMSTree.fjet1QGtagSub2_   = fjet->SubJet(1)->QGTag();
        fMitDMSTree.fjet1QGPtDSub2_   = fjet->SubJet(1)->QGPtD();
        fMitDMSTree.fjet1QGAxis1Sub2_ = fjet->SubJet(1)->QGAxis1();
        fMitDMSTree.fjet1QGAxis2Sub2_ = fjet->SubJet(1)->QGAxis2();
        fMitDMSTree.fjet1QGMultSub2_  = fjet->SubJet(1)->QGMult();
      }      

      // Angular info
      fMitDMSTree.fjet1metDphi_ = MathUtils::DeltaPhi(fjet->Phi(),(double)fMitDMSTree.metPhi_);
      fMitDMSTree.fjet1jet2Dphi_ = GetJetJetsDphi(fMitDMSTree.fjet1_);

    }// end filling of first fat jet

    // Increment cleaned fat jet counter
    fMitDMSTree.nfjets_++;

  }
     
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
    getGenLevelInfo(fMitDMSTree);

    Q         = fMCEventInfo->Scale();
    id1       = fMCEventInfo->Id1();
    x1        = fMCEventInfo->X1();
    pdf1      = fMCEventInfo->Pdf1();
    id2       = fMCEventInfo->Id2();
    x2        = fMCEventInfo->X2();
    pdf2      = fMCEventInfo->Pdf2();
    processId = fMCEventInfo->ProcessId();
  }

  fMitDMSTree.Q_ = Q;
  fMitDMSTree.id1_ = id1;
  fMitDMSTree.x1_ = x1;
  fMitDMSTree.pdf1_ = pdf1;
  fMitDMSTree.id2_ = id2;
  fMitDMSTree.x2_ = x2;
  fMitDMSTree.pdf2_ = pdf2;
  fMitDMSTree.processId_ = processId;

  // PILEUP RELATED
  if (! fIsData) {
    // loop over the pileup summary info and grab what you need
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing() ==  0)
        fMitDMSTree.npu_ = puinfo->GetPU_NumMean();
      if (puinfo->GetBunchCrossing() ==  1)
        fMitDMSTree.npuPlusOne_ = puinfo->GetPU_NumInteractions();
      if (puinfo->GetBunchCrossing() == -1)
        fMitDMSTree.npuMinusOne_ = puinfo->GetPU_NumInteractions();
    }
    fMitDMSTree.puweight_ = PUWeight(fMitDMSTree.npu_);
  }

  // Finally fill the tree
  fMitDMSTree.tree_->Fill();

  return;
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we request all
  // relevant information.

  if (! fIsData) {
    ReqBranch(fPileUpName,           fPileUp);
    ReqBranch(fMCEventInfoName,      fMCEventInfo);
    ReqEventObject(fMCParticlesName, fMCParticles,   true);
  }
  ReqEventObject(fPileUpDenName,     fPileUpDen,     true);
  ReqEventObject(fPVName,            fPV,            fPVFromBranch);
  ReqEventObject(fEvtSelDataName,    fEvtSelData,    true);
  ReqEventObject(fTriggerObjectsName,fTrigObj,       true);

  ReqEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch);
  ReqEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  ReqEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  ReqEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  ReqEventObject(fJetsName,          fJets,          fJetsFromBranch);
  ReqEventObject(fFatJetsName,       fFatJets,       fFatJetsFromBranch);
  ReqEventObject(fSubJetsName,       fSubJets,       fSubJetsFromBranch);
  ReqEventObject(fRawMetName,        fRawMet,        true);

  // Initialize the PU histrograms and weights
  // some useful definitions
  if (! fIsData) {
    TString dirFwk("AnaFwkMod");
    TString allEvts("hDAllEvents");
    // get input histo
    TFile *fif = new TFile(fPUInputFileName.Data());
    if (fif->IsOpen() == kFALSE) {
      printf(" WARNING -- missing input pile up file!\n");
    }  
    TDirectory *dirTmp = (TDirectory*) gROOT->FindObject(dirFwk.Data());
    if (dirTmp) {
      fif->cd(dirFwk.Data());
      fPUInput = (TH1D*) dirTmp->Get("hNPUTrue")->Clone();
      if (! fPUInput)
        printf(" WARNING -- no input framework file!\n");      
    }
    // get target histo
    TFile *ftf = new TFile(fPUTargetFileName.Data());
    if (ftf->IsOpen() == kFALSE) {
      printf(" WARNING -- missing target pile up file!\n");
    }  
    fPUTarget = (TH1D*) ftf->Get("pileup")->Clone();
    if (! fPUTarget)
      printf(" WARNING -- no target pile up histogram !\n");      
    // build pile up weight histo
    fPUInput->Rebin(10);
    fPUInput->Scale(1.0/fPUInput->GetSumOfWeights());
    fPUTarget->Scale(1.0/fPUTarget->GetSumOfWeights());
    fPUWeight = new TH1D((*fPUTarget) / (*fPUInput));
  }

  // Initialize DiJet MVA
  fDiJetMVA = new DiJetMVA();
  fDiJetMVA->Initialize(
          TString(getenv("CMSSW_BASE")+std::string("/src/MitMonoJet/Utils/data/vmva_weights/vtraining_lowpt_cen.root_BDTG.weights.xml")),
          TString(getenv("CMSSW_BASE")+std::string("/src/MitMonoJet/Utils/data/vmva_weights/vtraining_highpt_cen.root_BDTG.weights.xml")),
          std::string(getenv("CMSSW_BASE"))+std::string("/src/MitMonoJet/Utils/data/QGSystDatabase.txt"));          

  // Initialize Jet Uncertainties
  std::string jetCorrectorParams;
  std::string fatJetCorrectorParams;
  if (fIsData) {
    jetCorrectorParams = std::string(getenv("CMSSW_BASE"))+std::string("/src/MitPhysics/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt");
    fatJetCorrectorParams = std::string(getenv("CMSSW_BASE"))+std::string("/src/MitPhysics/data/FT53_V21A_AN6_Uncertainty_AK7PFchs.txt");
  }
  else {
    jetCorrectorParams = std::string(getenv("CMSSW_BASE"))+std::string("/src/MitPhysics/data/Summer13_V1_MC_Uncertainty_AK5PF.txt");
    fatJetCorrectorParams = std::string(getenv("CMSSW_BASE"))+std::string("/src/MitPhysics/data/FT53_V21A_AN6_Uncertainty_AK7PFchs.txt");
  }
  JetCorrectorParameters param(jetCorrectorParams);
  JetCorrectorParameters paramFat(fatJetCorrectorParams);
  fJetUncertainties = new JetCorrectionUncertainty(param);
  fFatJetUncertainties = new JetCorrectionUncertainty(paramFat);
  
  // Create Ntuple Tree
  fOutputFile = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  fMitDMSTree.CreateTree(0);
  fMitDMSTree.tree_->SetAutoSave(300e9);
  fMitDMSTree.tree_->SetDirectory(fOutputFile);
  AddOutput(fMitDMSTree.tree_);
}

//--------------------------------------------------------------------------------------------------
Float_t DMSTreeWriter::PUWeight(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!fPUWeight)
    return 1.0;
  
  return fPUWeight->GetBinContent(fPUWeight->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::getGenLevelInfo(MitDMSTree& tree)
{
  // Prepare 4-vector for invisible objects
  FourVector momInv;
  
  // Loop on all stable MC particles
  for (UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle *p = fMCParticles->At(i);
    if (p->Status()!=3)
      continue;

    // Check for Higgs(25),DM(100022),nuTau'(18) 
    if (p->AbsPdgId() == 25 || p->AbsPdgId() == 100022 || p->AbsPdgId() == 18)
      momInv += p->Mom();
    
    // Check if the particle is a Boson
    if (p->Is(MCParticle::kZ) || p->Is(MCParticle::kW)) {
      tree.genV_   = p->Mom();      
      tree.genVid_ = p->AbsPdgId();      

      // Check the daughters of the Boson
      if (p->NDaughters() > 0)
        tree.genVdaughterId_ = p->Daughter(0)->AbsPdgId();
        
      // Special case: top->Wb. Check if the mother is a top
      if (p->Mother()->Is(MCParticle::kTop)) {
        if (p->Mother()->PdgId()==6)      
          tree.topPt_    = p->Pt();
        else 
          tree.topBarPt_ = p->Pt();
      } // end top scope
      
    } // end boson scope 
    // Check if the particle is a high Pt photon
    else if (p->Is(MCParticle::kGamma) && p->Pt() > 100.) {
      tree.genV_   = p->Mom();      
      tree.genVid_ = p->AbsPdgId();            
    }
    else
      continue;
    
  } // end loop on MC Particles

  // Apply ttbar correction
  // reference is https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#MC_SFs_Reweighting
  if (tree.topPt_ > 0 && tree.topBarPt_ > 0) {
    float pt1 = TMath::Min((float)400., tree.topPt_);
    float pt2 = TMath::Min((float)400., tree.topBarPt_);
    float w1 = exp(0.156 - 0.00137*pt1);
    float w2 = exp(0.156 - 0.00137*pt2);
    tree.genweight_ = 1.001*sqrt(w1*w2); 
  }

  // Finish invisible object pt computation
  tree.genmet_ = momInv.Pt();
  tree.genmetPhi_ = momInv.Phi();

  return;
}
 
//--------------------------------------------------------------------------------------------------
Bool_t DMSTreeWriter::IsHLTMatched(LorentzVector& v,
                                   TriggerObject::ETriggerObject type,
                                   Float_t deltaR,
                                   Float_t minPt,
                                   Float_t maxEta)
{
  float minDr = 999.;
  
  // Loop on the selected trigger object collection 
  // and find a matched object
  for (UInt_t i=0; i<fTrigObj->GetEntries(); ++i) {
    const TriggerObject *to = fTrigObj->At(i);
    if (to->TriggerType() != type)
      continue; 
    // check trigger object pt,eta
    if (to->Pt() < minPt || fabs(to->Eta()) > maxEta)
      continue; 
    // compute the dR
    float thisDr = MathUtils::DeltaR(v, *to);
    if (thisDr < minDr) 
      minDr = thisDr;
    // check if user match condition is met
    if (minDr < deltaR)
      return true;
  } // end loop on trigger objects

  return false;
}


//--------------------------------------------------------------------------------------------------
Bool_t DMSTreeWriter::fjetIsCleaned(LorentzVector& v,
                                    Float_t deltaR)
{
  // Loop on tight muons and photons and discard 
  // jet if overlapping
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    // Discard loose muons
    if (fMuons->At(i)->ParticleId() == XsIsoParticle::EParticleId::eX)
      continue;
    float thisDr = MathUtils::DeltaR(v, fMuons->At(i)->Mom());
    if (thisDr < deltaR)
      return false;      
  }
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    // Discard soft photons
    if (fPhotons->At(i)->Pt() < 150.)
      continue;
    float thisDr = MathUtils::DeltaR(v, fPhotons->At(i)->Mom());
    if (thisDr < deltaR)
      return false;      
  }
 
  return true;      
}

//--------------------------------------------------------------------------------------------------
Int_t DMSTreeWriter::JetPartonMatch(LorentzVector& v,
                                    Float_t deltaR)
{
  float minDr = 999.;
  unsigned int pId = 0;
      
  // Loop on all stable MC particles and perform the matching
  for (UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle *p = fMCParticles->At(i);
    if (p->Status()!=3)
      continue;
    // compute the dR    
    float thisDr = MathUtils::DeltaR(v, *p);
    // standard check for any parton    
    if (thisDr < minDr) {
      minDr = thisDr;
      pId = p->AbsPdgId();
    }
    // special check: if quark inside matching cone and mother is W/Z
    // make extra checks and if true return the mother id
    if (thisDr < deltaR && p->AbsPdgId() < 6 
     && (p->Mother()->Is(MCParticle::kZ) || p->Mother()->Is(MCParticle::kW))) {
      // check that both the quarks are inside the matching cone
      float q1Dr = MathUtils::DeltaR(v, *(p->Mother()->Daughter(0)));
      float q2Dr = MathUtils::DeltaR(v, *(p->Mother()->Daughter(1)));
      if (q1Dr < deltaR && q2Dr < deltaR)
        return p->Mother()->AbsPdgId();
    }
  }

  // return matched parton only if minDr less that user dr limit
  if (minDr < deltaR)
    return pId;
  else 
    return 0;  
}

//--------------------------------------------------------------------------------------------------
Float_t DMSTreeWriter::GetJetJetsDphi(LorentzVector& v)
{
  // Loop on first two jets and consider the farthest one in phi  
  // which is not overlapping with the input jet (DR = 0.5)  
  UInt_t nMaxJets = 2;
  UInt_t nCleanJets = 0;
  float maxDphi = -10.;
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {
    const XlJet *jet = fJets->At(i);
    // Perform jet cleaning according to lepton counting
    if (fMitDMSTree.nlep_ > 0 
     && MathUtils::DeltaR(fMitDMSTree.lep1_,jet->Mom()) < 0.5)
      continue;
    if (fMitDMSTree.nlep_ > 1 
     && MathUtils::DeltaR(fMitDMSTree.lep2_,jet->Mom()) < 0.5)
      continue;
    if (nCleanJets >= nMaxJets) 
      continue;
    nCleanJets++;
    // Discard ovelapping jets
    if (MathUtils::DeltaR(v, jet->Mom()) < 0.5)
      continue;
    float thisDphi = MathUtils::DeltaPhi(v.Phi(),jet->Phi());
    if (thisDphi > maxDphi)
      maxDphi = thisDphi;
  }
 
  return maxDphi;      
}

//--------------------------------------------------------------------------------------------------
Float_t DMSTreeWriter::GetFatJetBtag(LorentzVector& v,
                                     Float_t deltaR)
{
  float minDr = 999.;
    
  // Loop on standard jet collection, find the best matched jet 
  // and get the btagging
  for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
    const Jet *jet = fJets->At(i);
    // compute the dR
    float thisDr = MathUtils::DeltaR(v, *jet);
    if (thisDr < minDr) 
      minDr = thisDr;
    // if the fat jet is matched to standard jet return the btagging
    // of the standard jet
    if (minDr < deltaR)
      return jet->CombinedSecondaryVertexBJetTagsDisc();  
  }
      
  return -1;
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::CorrectMet(const float met, const float metPhi, const LorentzVector& l1, const LorentzVector& l2,
                               float& newMet, float& newMetPhi)
{
  // inputs:  met, metPhi, l1, l2
  // outputs: newMet, newMetPhi
  float newMetX;
  float newMetY;
  newMetX = met*TMath::Cos(metPhi) + l1.Px() + l2.Px();
  newMetY = met*TMath::Sin(metPhi) + l1.Py() + l2.Py();
  
  newMet    = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);

  return;
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::CorrectMet(const float met, const float metPhi, const LorentzVector& l1,
                               float& newMet, float& newMetPhi, Bool_t applyFprintCorrection)
{
  // inputs:  met, metPhi, l1
  // outputs: newMet, newMetPhi
  float newMetX;
  float newMetY;
  float fPrintCorr = 0.;
  if (applyFprintCorrection)
    fPrintCorr = GetFootprint(l1.Pt(),l1.Eta());  
  
  newMetX = met*TMath::Cos(metPhi) + (l1.Pt()+fPrintCorr)*TMath::Cos(l1.Phi());
  newMetY = met*TMath::Sin(metPhi) + (l1.Pt()+fPrintCorr)*TMath::Sin(l1.Phi());
    
  newMet    = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);

  return;
}

//--------------------------------------------------------------------------------------------------
float DMSTreeWriter::GetMt(const LorentzVector& l1, const float met, const float metPhi)
{
  // inputs:  l1, met, metPhi
  // outputs: transverse mass
  double lepPhi = l1.Phi();
  return TMath::Sqrt(2*met*l1.Pt()*(1-TMath::Cos(MathUtils::DeltaPhi(lepPhi,(double)metPhi))));
}

//--------------------------------------------------------------------------------------------------
float DMSTreeWriter::GetFootprint(const float pt, const float eta)
{
  // inputs:  pt, eta
  // outputs: footprint pt value
  if (fIsData) {
    if (fabs(eta) < 1.5) 
      return (0.4633 + TMath::Min(pt,(float)130.)*0.0055);
    else
      return (0.3951 + TMath::Min(pt,(float)140.)*0.0125); //140 is NOT a typo!
  }
  else {
    if (fabs(eta) < 1.5) 
      return (0.5090 + TMath::Min(pt,(float)130.)*0.0075);
    else
      return (0.2969 + TMath::Min(pt,(float)130.)*0.0135);
  }
    
}
