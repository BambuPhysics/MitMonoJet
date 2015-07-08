  // $Id $

#include <iostream>
#include <sstream>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"

#include "MitMonoJet/SelMods/interface/BoostedVAnalysisMod.h"


using namespace mithep;
ClassImp(mithep::BoostedVAnalysisMod)

//--------------------------------------------------------------------------------------------------
  BoostedVAnalysisMod::BoostedVAnalysisMod(const char *name, const char *title) :
    BaseMod(name,title),
    // define all the Branches to load
    fMetBranchName         ("PFMet"),
    fFatJetsName           ("AKt7PFJets"),
    fJetsName              (Names::gkPFJetBrn),
    fElectronsName         (Names::gkElectronBrn),
    fMuonsName             (Names::gkMuonBrn),
    fPhotonsName           (Names::gkPhotonBrn),
    fTriggerObjectsName    ("HltObjsMonoJet"),
    fMetFromBranch         (kTRUE),
    fFatJetsFromBranch     (kTRUE),
    fJetsFromBranch        (kTRUE),
    fElectronsFromBranch   (kTRUE),
    fMuonsFromBranch       (kTRUE),
    fPhotonsFromBranch     (kTRUE),
    // define active preselection regions
    fApplyResolvedPresel   (kTRUE),
    fApplyTopPresel        (kTRUE),
    fApplyWlepPresel       (kTRUE),
    fApplyZmmPresel        (kTRUE),
    fApplyZeePresel        (kTRUE),
    fApplyMetPresel        (kTRUE), 
    fApplyVbfPresel        (kTRUE),
    fApplyGjetPresel       (kTRUE),
    fApplyFatJetPresel     (kTRUE),
    // define other flags
    fFillAndPublishPresel  (kTRUE),
    fSkipEvents            (kTRUE),
    // collections
    fMet                   (0),
    fJets                  (0),
    fElectrons             (0),
    fMuons                 (0),
    fPhotons               (0),
    fEvtSelData            (0),    
    fXlEvtSelData          (0),
    // cuts
    fMinResolvedMass       (60),
    fMaxResolvedMass       (110),
    fMinFatJetPt           (200),
    fMinTagJetPt           (100),
    fMinVbfJetPt           (40),
    fMinMet                (200),
    fMinVbfMass            (800),
    fMinVbfMet             (110),
    fMinPhotonPt           (160),
    // counters
    fAll                   (0),
    fPass                  (0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
BoostedVAnalysisMod::~BoostedVAnalysisMod()
{
  // Destructor
  delete fXlEvtSelData;
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::SlaveBegin()
{
  // Load Branches
  // ReqEventObject(fMetBranchName, fMet,       fMetFromBranch);
  // ReqEventObject(fJetsName,      fJets,      fJetsFromBranch);
  // ReqEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  // ReqEventObject(fMuonsName,     fMuons,     fMuonsFromBranch);
  // ReqEventObject(fPhotonsName,   fPhotons,   fPhotonsFromBranch);
  
  // If requested by the user prepare the collections to store the
  // preselection word
  if (fFillAndPublishPresel) {
    ReqEventObject(Names::gkEvtSelDataBrn, fEvtSelData, true);
    // Create the new output collection
    fXlEvtSelData = new XlEvtSelData("XlEvtSelData");
    PublishObj(fXlEvtSelData);
  }

}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Process()
{
  // LoadEventObject(fMetBranchName, fMet,       fMetFromBranch);
  // LoadEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  // LoadEventObject(fMuonsName,     fMuons,     fMuonsFromBranch);
  // LoadEventObject(fPhotonsName,   fPhotons,   fPhotonsFromBranch);
  
  fElectrons = GetObject<ElectronOArr>(fElectronsName);
  fMet = GetObject<PFMetCol>(fMetBranchName);
  fMuons = GetObject<MuonOArr>(fMuonsName);
  fPhotons = GetObject<PhotonOArr>(fPhotonsName);

  if (fFillAndPublishPresel) {
    LoadEventObject(Names::gkEvtSelDataBrn, fEvtSelData, true);
  }

  fJets    = GetObject<JetOArr>(fJetsName);
  if (fApplyFatJetPresel)
    fFatJets = GetObject<JetOArr>(fFatJetsName);

  // Setup HLT bits
  Bool_t passSingleMuHLT = kTRUE;
  Bool_t passMonoJetHLT = kTRUE;
  Bool_t passVbfHLT = kTRUE;
  Bool_t passGjetHLT = kTRUE;
  Bool_t passDiEleHLT = kTRUE;
  fTrigObj = GetHLTObjects(fTriggerObjectsName);
  if (! fTrigObj)
   // printf("BoostedVAnalysisMod::TriggerObjectCol not found\n");
   ;
  else {
    // Loop through the stored trigger objects and find corresponding trigger name
    for (UInt_t i=0;i<fTrigObj->GetEntries();++i) {
      const TriggerObject *to = fTrigObj->At(i);
      TString trName = to->TrigName();
      // Default Single muon
      if (trName.Contains("HLT_IsoMu24") ||
          trName.Contains("HLT_IsoMu17") ||
          trName.Contains("HLT_IsoMu15"))
        passSingleMuHLT = kTRUE;
      // Default MonoJet
      if (trName.Contains("MonoCentralPFJet80_PFMETnoMu") ||
          trName.Contains("HLT_MET120_HBHENoiseCleaned_v"))
        passMonoJetHLT = kTRUE;
      // Default VBF
      if (trName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"))
        passVbfHLT = kTRUE;
      // Default Photon+jets
      if (trName.Contains("HLT_Photon135") ||
          trName.Contains("HLT_Photon150"))
        passGjetHLT = kTRUE;
      // Default DiEle
      if (trName.Contains("HLT_Ele17"))
        passDiEleHLT = kTRUE;
    }
  }
   
  // Initialize selection flags
  Bool_t passResolvedPresel = kFALSE;
  Bool_t passTopPresel = kFALSE;
  Bool_t passWlepPresel = kFALSE;
  Bool_t passZmmPresel = kFALSE;
  Bool_t passZeePresel = kFALSE;
  Bool_t passMetPresel = kFALSE;
  Bool_t passVbfPresel = kFALSE;
  Bool_t passGjetPresel = kFALSE;

  // Increment all events counter
  fAll++;
  
  // Discard events with no jets
  if (fJets->GetEntries() < 1 && fSkipEvents) {
    this->SkipEvent(); 
    return;
  }

  // Determine if event passes resolved preselection (di-jet, di-jet mass)
  if (fApplyResolvedPresel && fJets->GetEntries() > 1) {
    int nGoodJetPairs = 0;

    // Jets, check btagging for vetoing. Break loop as soon as one pair is found
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jetOne = fJets->At(i);
      // Pt and eta cuts
      if (jetOne->Pt() < 30. || fabs(jetOne->Eta()) > 2.5)
        continue;
      // Check btagging
      if (jetOne->CombinedSecondaryVertexBJetTagsDisc() > 0.679)
        continue;
      // di-jet mass cut                
      for (UInt_t j = i+1; j < fJets->GetEntries(); ++j) {
        const Jet *jetTwo = fJets->At(j);
        // Pt and eta cuts
        if (jetTwo->Pt() < 30. || fabs(jetTwo->Eta()) > 2.5)
          continue;
        // Check btagging
        if (jetTwo->CombinedSecondaryVertexBJetTagsDisc() > 0.679)
          continue;
        // Check mass
        if ((jetOne->Mom() + jetTwo->Mom()).M() > fMinResolvedMass
          &&(jetOne->Mom() + jetTwo->Mom()).M() < fMaxResolvedMass) {
          nGoodJetPairs++;
          break;
        } // end loop on second jet
      }
      if (nGoodJetPairs > 0)
        break;        
    }
    
    // Decision    
    if (nGoodJetPairs > 0)
      passResolvedPresel = kTRUE;
  }

  // Determine if the event passes any of the preselection cuts
  if (fApplyTopPresel) {
    // Top Preselection: require boosted jet + 2 b-jets + lepton
    // Following TOP-12-042 PAS selections
    // also require standard tag jet selection for trigger matching
    int nGoodFatJets = 0;
    int nGoodTagJets = 0;
    int nGoodBJets = 0;
    int nGoodLeptons = 0;

    // FatJets: if fatJets not available adjust preselection params
    if (!fApplyFatJetPresel)
      nGoodFatJets = 1;
    if (fApplyFatJetPresel && fFatJets->GetEntries() > 0) { 
      for (UInt_t i = 0; i < fFatJets->GetEntries(); ++i) {
        const Jet *jet = fFatJets->At(i);
        // Pt and eta cuts
        if (jet->Pt() < fMinFatJetPt || fabs(jet->Eta()) > 2.5)
          nGoodFatJets++;
      }
    }

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < 30. || fabs(jet->Eta()) > 2.5)
        continue;
      // Check high pt jet
      if (jet->Pt() > fMinTagJetPt)
        nGoodTagJets++;
      // Check b-tagging
      if (jet->CombinedSecondaryVertexBJetTagsDisc() > 0.244)
        nGoodBJets++;
    }

    // Leptons
    if (fElectrons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fElectrons->GetEntries(); ++i) {
        const Electron *ele = fElectrons->At(i);
        // Pt and eta cuts
        if (ele->Pt() < 30. || fabs(ele->Eta()) > 2.5 || (fabs(ele->Eta()) > 1.442 && fabs(ele->Eta()) > 1.566))
          continue;
        nGoodLeptons++;
      }
    }
    if (fMuons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fMuons->GetEntries(); ++i) {
        const Muon *mu = fMuons->At(i);
        // Pt and eta cuts
        if (mu->Pt() < 30. || fabs(mu->Eta()) > 2.1)
          continue;
        nGoodLeptons++;
      }
    }
    if (passSingleMuHLT && nGoodTagJets > 0 && nGoodBJets > 1 && nGoodFatJets > 0 && nGoodLeptons > 0)
      passTopPresel = kTRUE;
  }

  if (fApplyWlepPresel) {
    // W Preselection: require high pt jet/resolved pair + muon
    int nGoodTagJets = 0;
    int nGoodLeptons = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinTagJetPt || fabs(jet->Eta()) > 2.5)
        continue;
      nGoodTagJets++;
    }

    // Leptons
    if (fMuons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fMuons->GetEntries(); ++i) {
        const Particle *lep = fMuons->At(i);
        // Pt cut
        if (lep->Pt() < 15.)
          continue;
        nGoodLeptons++;
      }
    }
            
    if (passMonoJetHLT && nGoodTagJets > 0 && nGoodLeptons > 0)
      passWlepPresel = kTRUE;
    if (passMonoJetHLT && passResolvedPresel && nGoodLeptons > 0)
      passWlepPresel = kTRUE;
  }

  if (fApplyZmmPresel || fApplyZeePresel) {
    // Zmm Preselection: require high pt jet/resolved pair + di-leptons
    int nGoodTagJets = 0;
    int nGoodMuons = 0;
    int nGoodEles = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinTagJetPt || fabs(jet->Eta()) > 2.5)
        continue;
      nGoodTagJets++;
    }

    // Leptons
    if (fMuons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fMuons->GetEntries(); ++i) {
        const Particle *lep = fMuons->At(i);
        // Pt cut
        if (lep->Pt() < 10.)
          continue;
        nGoodMuons++;
      }
    }
    if (fElectrons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fElectrons->GetEntries(); ++i) {
        const Particle *lep = fElectrons->At(i);
        // Pt cut
        if (lep->Pt() < 10.)
          continue;
        nGoodEles++;
      }
    }
            
    if (passMonoJetHLT && nGoodTagJets > 0 && nGoodMuons > 1)
      passZmmPresel = kTRUE;
    if (passMonoJetHLT && passResolvedPresel && nGoodMuons > 1)
      passZmmPresel = kTRUE;
    if (passDiEleHLT && nGoodTagJets > 0 && nGoodEles > 1)
      passZeePresel = kTRUE;
    if (passDiEleHLT && passResolvedPresel && nGoodEles > 1)
      passZeePresel = kTRUE;
  }
            
  if (fApplyMetPresel) {
    // Met Preselection: require boosted jet/resolved pair + MET
    int nGoodTagJets = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinTagJetPt || fabs(jet->Eta()) > 2.5)
        continue;
      nGoodTagJets++;
    }
        
    if (passMonoJetHLT && nGoodTagJets > 0 && fMet->At(0)->Pt() > fMinMet)
      passMetPresel = kTRUE;
    if (passMonoJetHLT && passResolvedPresel && fMet->At(0)->Pt() > fMinMet)
      passMetPresel = kTRUE;
  }

  if (fApplyVbfPresel) {
    // Vbf Preselection: require VBF jets + METnoMu
    int nGoodVbfPairs = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jetOne = fJets->At(i);
      // Pt and eta cuts
      if (jetOne->Pt() < fMinVbfJetPt)
        continue;
      // Jet mulitplicity cut (at least need second jet)
      if (fJets->GetEntries() < 2) 
        continue;
      // di-jet mass cut                
      for (UInt_t j = i+1; j < fJets->GetEntries(); ++j) {
        const Jet *jetTwo = fJets->At(j);
        if (jetTwo->Pt() < fMinVbfJetPt) 
          continue;
        if ((jetOne->Mom() + jetTwo->Mom()).M() > fMinVbfMass)
          nGoodVbfPairs++;
      }
                
    }
    
    // METnoMu
    float METnoMuPt, METnoMuPhi;
    if (fMuons->GetEntries() > 1) {
      CorrectMet(fMet->At(0)->Pt(), fMet->At(0)->Phi(),fMuons->At(0),fMuons->At(1),
                 METnoMuPt,METnoMuPhi);
    }
    else if (fMuons->GetEntries() == 1) {
      CorrectMet(fMet->At(0)->Pt(), fMet->At(0)->Phi(),fMuons->At(0),0,
                 METnoMuPt,METnoMuPhi);
    }
    else 
      METnoMuPt = fMet->At(0)->Pt();                      
        
    if (passVbfHLT && nGoodVbfPairs > 0 && METnoMuPt > fMinVbfMet)
      passVbfPresel = kTRUE;
  }

  if (fApplyGjetPresel) {
    // G+jets Preselection: require tag jet/resolved pair + photon
    int nGoodTagJets = 0;
    int nGoodPhotons = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinTagJetPt || fabs(jet->Eta()) > 2.5)
        continue;
      nGoodTagJets++;
    }

    // Photons
    if (fPhotons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fPhotons->GetEntries(); ++i) {
        const Photon *pho = fPhotons->At(i);
        // Pt and eta cuts
        if (pho->Pt() < fMinPhotonPt || fabs(pho->Eta()) > 2.5)
          continue;
        nGoodPhotons++;
      }
    }
    if (passGjetHLT && nGoodTagJets > 0 && nGoodPhotons > 0)
      passGjetPresel = kTRUE;
    if (passGjetHLT && passResolvedPresel && nGoodPhotons > 0)
      passGjetPresel = kTRUE;
  }

  // Skip event if it does not pass any preselection
  if (!passTopPresel && !passWlepPresel && !passZmmPresel && !passZeePresel 
   && !passMetPresel && !passVbfPresel  && !passGjetPresel
   && fSkipEvents) {
    this->SkipEvent(); 
    return;
  }

  // Store the preselection word in the extended event selection object
  if (fFillAndPublishPresel) {
    int HLTWord          = GetHLTWord (
                           passSingleMuHLT, 
                           passMonoJetHLT,
                           passVbfHLT,
                           passGjetHLT,
                           passDiEleHLT); 
    int preselectionWord = GetPreselWord (
                           passTopPresel, 
                           passWlepPresel,
                           passZmmPresel,
                           passZeePresel,
                           passMetPresel, 
                           passVbfPresel, 
                           passGjetPresel,
                           passResolvedPresel); 
    fXlEvtSelData->SetFiltersWord(fEvtSelData->metFiltersWord());
    fXlEvtSelData->SetHLTWord(HLTWord);
    fXlEvtSelData->SetPreselWord(preselectionWord);
  }
   
  // Increment passed events counter
  fPass++;
  return;
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::SlaveTerminate()
{

  // Terminate preselection and print out module selection efficiency
  Double_t frac =  100.*fPass/fAll;
  Info("SlaveTerminate", "Selected %.2f%% events (%lld out of %lld)", 
       frac, fPass, fAll);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::CorrectMet(const float met, const float metPhi,
                                     const Particle *l1, const Particle *l2,
                                     float &newMet, float &newMetPhi)
{
  // inputs: met, metPhi, l1, [ l2 only used if pointer is non-zero ]
  // outputs: newMet, newMetPhi

  float newMetX;
  float newMetY;
  if (l2) { // these are doubles on the right: need full calculation in one line
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px() + l2->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py() + l2->Mom().Py();
  }
  else {
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py();
  }
  
  newMet = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);
}

//--------------------------------------------------------------------------------------------------
int  BoostedVAnalysisMod::GetHLTWord( 
                          bool passSingleMuHLT,
                          bool passMonoJetHLT,
                          bool passVbfHLT,
                          bool passGjetHLT,
                          bool passDiEleHLT)
{  
  // This function creates the word containing the HLT bit decisions.
  // The bit ordering follows the order of the parameters passed
  // to this function. 
  
  //Initialize the word
  int theWord = 0;
  //Initialize the vector of bits
  std::vector<int> theBits;
  theBits.push_back((int) passSingleMuHLT);
  theBits.push_back((int) passMonoJetHLT);
  theBits.push_back((int) passVbfHLT);
  theBits.push_back((int) passGjetHLT);
  theBits.push_back((int) passDiEleHLT);
  //Create the word
  for (unsigned int iBit = 0; iBit < theBits.size(); iBit++)
    theWord |= theBits[iBit] << iBit;
  
  return theWord;
}

//--------------------------------------------------------------------------------------------------
int  BoostedVAnalysisMod::GetPreselWord( 
                          bool passTopPresel, 
                          bool passWlepPresel,
                          bool passZmmPresel,
                          bool passZeePresel,
                          bool passMetPresel, 
                          bool passVbfPresel,
                          bool passGjetPresel,
                          bool passResolvedPresel)
{  
  // This function creates the word containing the bit decisions.
  // The bit ordering follows the order of the parameters passed
  // to this function. 
  
  //Initialize the word
  int theWord = 0;
  //Initialize the vector of bits
  std::vector<int> theBits;
  theBits.push_back((int) passTopPresel);
  theBits.push_back((int) passWlepPresel);
  theBits.push_back((int) passZmmPresel);
  theBits.push_back((int) passZeePresel);
  theBits.push_back((int) passMetPresel);
  theBits.push_back((int) passVbfPresel);
  theBits.push_back((int) passGjetPresel);
  theBits.push_back((int) passResolvedPresel);
  //Create the word
  for (unsigned int iBit = 0; iBit < theBits.size(); iBit++)
    theWord |= theBits[iBit] << iBit;
  
  return theWord;
}
  
