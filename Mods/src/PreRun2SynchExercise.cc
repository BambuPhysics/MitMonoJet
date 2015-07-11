#include "MitMonoJet/Mods/interface/PreRun2SynchExercise.h"

#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFJet.h"

#include "TVector2.h"

ClassImp(mithep::PreRun2SynchExercise)

void
mithep::PreRun2SynchExercise::Process()
{
  auto* inVertices = GetObject<mithep::VertexCol>(fVerticesName);
  auto* inMets = GetObject<mithep::MetCol>(fMetName);
  auto* inJets = GetObject<mithep::JetCol>(fJetsName);
  auto* inElectrons = GetObject<mithep::ElectronCol>(fElectronsName);
  auto* inMuons = GetObject<mithep::MuonCol>(fMuonsName);
  auto* inTaus = GetObject<mithep::PFTauCol>(fTausName);
  auto* inPhotons = GetObject<mithep::PhotonCol>(fPhotonsName);

  unsigned nTaus = 0;
  for (unsigned iT = 0; iT != inTaus->GetEntries(); ++iT) {
    if (inTaus->At(iT)->PFTauDiscriminator(mithep::PFTau::kDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits) < 5.)
      ++nTaus;
  }

  fSynchPass[kGoodVertex] = inVertices->GetEntries() != 0;
  fSynchPass[kMET200] = inMets->At(0)->Pt() > 200.;
  fSynchPass[kNAK4Jets] = inJets->GetEntries() < 3;
  fSynchPass[kEMuVeto] = inElectrons->GetEntries() + inMuons->GetEntries() == 0;
  fSynchPass[kTauVeto] = nTaus == 0;
  fSynchPass[kPhotonVeto] = inPhotons->GetEntries() == 0;

  fSynchPass[kCleanJets] = false;
  fSynchPass[kLeadingJet110] = false;
  fSynchPass[kDeltaPhiJ1J2] = false;

  if (inJets->GetEntries() > 0) {
    PFJet const* leadJet = static_cast<PFJet const*>(inJets->At(0));

    double rawE = leadJet->RawMom().E();

    if (leadJet->ChargedHadronEnergy() / rawE > 0.2 &&
        leadJet->NeutralHadronEnergy() / rawE < 0.7 &&
        leadJet->NeutralEmEnergy() / rawE < 0.7) {
      fSynchPass[kCleanJets] = true;
      if (leadJet->Pt() > 110.)
        fSynchPass[kLeadingJet110] = true;
    }

    fSynchPass[kDeltaPhiJ1J2] = true;

    if (inJets->GetEntries() > 1) {
      PFJet const* trailJet = static_cast<PFJet const*>(inJets->At(1));

      double rawETrail = trailJet->RawMom().E();

      if (trailJet->NeutralHadronEnergy() / rawETrail < 0.7 &&
          trailJet->NeutralEmEnergy() / rawETrail < 0.9) {
        if (std::abs(TVector2::Phi_mpi_pi(leadJet->Phi() - trailJet->Phi())) > 2.5)
          fSynchPass[kDeltaPhiJ1J2] = false;
      }
      else {
        fSynchPass[kCleanJets] = false;
        fSynchPass[kDeltaPhiJ1J2] = false;
      }
    }
  }

  fTree->Fill();
}

void
mithep::PreRun2SynchExercise::SlaveBegin()
{
  fTree = new TTree("synchTree", "Synch checkpoints");
  fTree->Branch("GoodVertex", fSynchPass + kGoodVertex, "pass/O");
  fTree->Branch("LeadingJet110", fSynchPass + kLeadingJet110, "pass/O");
  fTree->Branch("DeltaPhiJ1J2", fSynchPass + kDeltaPhiJ1J2, "pass/O"); 
  fTree->Branch("MET200", fSynchPass + kMET200, "pass/O");       
  fTree->Branch("NAK4Jets", fSynchPass + kNAK4Jets, "pass/O");     
  fTree->Branch("EMuVeto", fSynchPass + kEMuVeto, "pass/O");      
  fTree->Branch("TauVeto", fSynchPass + kTauVeto, "pass/O");      
  fTree->Branch("PhotonVeto", fSynchPass + kPhotonVeto, "pass/O");   

  AddOutput(fTree);
}
