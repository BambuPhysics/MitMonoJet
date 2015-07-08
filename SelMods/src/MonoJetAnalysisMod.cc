// $Id $

#include <iostream>
#include <sstream>
#include <limits>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"

#include "MitMonoJet/SelMods/interface/MonoJetAnalysisMod.h"


using namespace mithep;
ClassImp(mithep::MonoJetAnalysisMod)

//--------------------------------------------------------------------------------------------------
  MonoJetAnalysisMod::MonoJetAnalysisMod(const char *name, const char *title) :
    BaseMod(name,title),
    //define all the Branches to load
    fCategoriesName        ("MonoJetCategories"),
    fMetBranchName         ("PFMet"),
    fJetsName              (Names::gkPFJetBrn),
    fElectronsName         (Names::gkElectronBrn),
    fMuonsName             (Names::gkMuonBrn),
    fTausName              (Names::gkPFTauBrn),
    fLeptonsName           (ModNames::gkMergedLeptonsName),
    fPFCandidatesName      (Names::gkPFCandidatesBrn),
    fMetFromBranch         (kTRUE),
    fJetsFromBranch        (kTRUE),
    fPFCandidatesFromBranch(kTRUE),
    fElectronsFromBranch   (kTRUE),
    fMuonsFromBranch       (kTRUE),
    fTausFromBranch        (kTRUE),
    // collections
    fMet                   (0),
    fJets                  (0),
    fElectrons             (0),
    fMuons                 (0),
    fPFTaus                (0),
    fPFCandidates          (0),
    // event category flags
    fCategories            (0),
    // counters
    fNEventsSelected       (0)
{
  // cuts
  UInt_t ibig = UInt_t(-1);
  Double_t dbig = std::numeric_limits<double>::max();
  std::fill_n(fMinNumLeptons, nCat, ibig);
  std::fill_n(fMaxNumLeptons, nCat, 0);
  std::fill_n(fMinNumTaus, nCat, ibig);
  std::fill_n(fMaxNumTaus, nCat, 0);
  std::fill_n(fMinNumJets, nCat, ibig);
  std::fill_n(fMaxNumJets, nCat, 0);
  std::fill_n(fMinNumGenNeutrinos, nCat, ibig);
  std::fill_n(fMinJetEt, nCat, dbig);
  std::fill_n(fMaxJetEta, nCat, 0.);
  std::fill_n(fMinMetEt, nCat, dbig);
  std::fill_n(fMinEmulMetEt, nCat, dbig);
  std::fill_n(fMinChargedHadronFrac, nCat, dbig);
  std::fill_n(fMaxNeutralHadronFrac, nCat, 0.);
  std::fill_n(fMaxNeutralEmFrac, nCat, 0.);
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::SlaveBegin()
{
  // Selection Histograms
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    // Histograms after preselection
    AddTH1(fMetEt[iCat], TString::Format("hMetEt%d", iCat), ";MetEt;Number of Events" ,400,0.,400.);
    AddTH1(fJetEt[iCat], TString::Format("hJetEt%d", iCat), ";JetEt;Number of Events" ,400,0.,400.);
    AddTH1(fJetEta[iCat],TString::Format("hJetEta%d", iCat),";JetEta;Number of Events",50,-5,5);
  }

  fCategories = new mithep::Array<mithep::TriggerMask>(1, fCategoriesName);
  fCategories->AddNew();
  fCategories->SetMustDeleteBit();
  PublishObj(fCategories);
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Process()
{
  fMet = GetObject<MetCol>(fMetBranchName);
  fElectrons = GetObject<ElectronCol>(fElectronsName);
  fMuons = GetObject<MuonCol>(fMuonsName);
  fPFTaus = GetObject<PFTauCol>(fTausName);
  fPFCandidates = GetObject<PFCandidateCol>(fPFCandidatesName);

  fJets = GetObject<JetOArr>(fJetsName);

  ParticleOArr   *leptons      = GetObject<ParticleOArr>(fLeptonsName);
  MCParticleOArr *genNeutrinos = GetObject<MCParticleOArr>(ModNames::gkMCNeutrinosName);

  // Define Cuts
  const int nCuts = 7;
  bool passCut[nCat][nCuts];
  for(unsigned iCat = 0; iCat != nCat; ++iCat)
    std::fill_n(passCut[iCat], nCuts, false);

  // Discard events with no identified required objects (jets, leptons, taus)
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    if (fJets->GetEntries() >= fMinNumJets[iCat] && fJets->GetEntries() <= fMaxNumJets[iCat])
      passCut[iCat][0] = true;
    if (leptons->GetEntries() >= fMinNumLeptons[iCat] && leptons->GetEntries() <= fMaxNumLeptons[iCat])
      passCut[iCat][1] = true;
    if (fPFTaus->GetEntries() >= fMinNumTaus[iCat] && fPFTaus->GetEntries() <= fMaxNumTaus[iCat])
      passCut[iCat][2] = true;
  }

  std::vector<int> theGoodJets[nCat];

  // Discard events with soft jets
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    int nHardJet = 0;
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      if (jet->Et() <= fMinJetEt[iCat])
        continue;
      if (TMath::Abs(jet->Eta()) >= fMaxJetEta[iCat])
        continue;

      nHardJet++;
      theGoodJets[iCat].push_back((int) i);
    }
    if (nHardJet > 0)
      passCut[iCat][3] = true;
  }

  // Discard events with soft met
  const Met *stdMet = fMet->At(0);
  mithep::FourVectorM emulMet = stdMet->Mom();
  for (unsigned iL = 0; iL != leptons->GetEntries(); ++iL)
    emulMet += leptons->At(iL)->Mom();
  
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    if (stdMet->Pt() > fMinMetEt[iCat])
      passCut[iCat][4] = true;

    if (emulMet.Pt() > fMinEmulMetEt[iCat])
      passCut[iCat][5] = true;
  }

  UInt_t nuInAcceptance = 0;
  for (UInt_t i = 0; i < genNeutrinos->GetEntries(); ++i)
    if (fabs(genNeutrinos->At(i)->Eta())<2.1 and genNeutrinos->At(i)->Pt()>10) nuInAcceptance++;

  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    if(nuInAcceptance >= fMinNumGenNeutrinos[iCat])
      passCut[iCat][6] = true;
  }

  mithep::BitMask1024 catBits;
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    bool passAllCuts = true;
    for (int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[iCat][c];
    if (passAllCuts) {
      fNEventsSelected++;
      // Make Preselection Histograms
      for (int iGoodJet=0; iGoodJet<(int)theGoodJets[iCat].size(); iGoodJet++)
        fJetEt[iCat]->Fill(fJets->At(theGoodJets[iCat][iGoodJet])->Et());
      fMetEt[iCat]->Fill(fMet->At(0)->Et());
      fJetEta[iCat]->Fill(fJets->At(0)->Eta());

      catBits.SetBit(iCat);
    }
  }

  if(catBits.NBitsSet() == 0){
    this->SkipEvent(); // skip the event if does not passes the cuts
  }

  fCategories->At(0)->SetBits(catBits);
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::SlaveTerminate()
{
  std::cout << "selected events on MonoJetAnalysisMod: " << fNEventsSelected << std::endl;

  fCategories->Delete();
}
//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Terminate()
{
}
