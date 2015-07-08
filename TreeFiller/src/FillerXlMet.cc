#include <TSystem.h> 

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlMet.h"
#include "MitAna/DataTree/interface/XlMet.h"

using namespace mithep;

ClassImp(mithep::FillerXlMet)

//--------------------------------------------------------------------------------------------------
FillerXlMet::FillerXlMet(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fJetsName (Names::gkPFJetBrn),
  fMuonsName (Names::gkMuonBrn),
  fElectronsName (Names::gkElectronBrn),
  fTausName (Names::gkPFTauBrn),
  fPhotonsName (Names::gkPhotonBrn),
  fPFCandidatesName (Names::gkPFCandidatesBrn),
  fPVName (Names::gkPVBeamSpotBrn),
  fPileUpDenName (Names::gkPileupEnergyDensityBrn),
  fRawMetName ("PFMet"),
  fXlMetName ("XlMetMVA"),
  fJetsFromBranch (kTRUE),
  fMuonsFromBranch (kTRUE),
  fElectronsFromBranch (kTRUE),
  fTausFromBranch (kTRUE),
  fPhotonsFromBranch (kTRUE),
  fPFCandidatesFromBranch (kTRUE),
  fPVFromBranch (kTRUE),
  fPublishOutput (kTRUE),
  fXlMet (0)
{
  // Constructor.
}

FillerXlMet::~FillerXlMet()
{
  // Destructor
  delete fXlMet;
}

//--------------------------------------------------------------------------------------------------
void FillerXlMet::Process()
{
  // make sure the out collections are empty before starting
  fXlMet->Delete();  
  
  // Load the branches we want to work with
  fJets = GetObject<JetOArr>(fJetsName);
  fMuons = GetObject<MuonOArr>(fMuonsName);
  fElectrons = GetObject<ElectronOArr>(fElectronsName);
  fPFTaus = GetObject<PFTauOArr>(fTausName);
  fPhotons = GetObject<PhotonOArr>(fPhotonsName);
  fPFCandidates = GetObject<PFCandidateCol>(fPFCandidatesName);
  fPV = GetObject<VertexCol>(fPVName);
  fRawMet = GetObject<PFMetCol>(fRawMetName);
  fPileUpDen = GetObject<PileupEnergyDensityCol>(fPileUpDenName);
  
  // LoadEventObject(fJetsName, fJets, fJetsFromBranch);
  // LoadEventObject(fMuonsName, fMuons, fMuonsFromBranch);
  // LoadEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  // LoadEventObject(fTausName, fPFTaus, fTausFromBranch);
  // LoadEventObject(fPhotonsName, fPhotons, fPhotonsFromBranch);
  // LoadEventObject(fPFCandidatesName,  fPFCandidates,  fPFCandidatesFromBranch);
  // LoadEventObject(fPVName, fPV, fPVFromBranch);
  // LoadEventObject(fRawMetName, fRawMet, true);
  // LoadEventObject(fPileUpDenName, fPileUpDen, true);

  // Convert the input collection into PFJets
  PFJetOArr *fPFJets = new PFJetOArr;
  for(UInt_t i=0; i<fJets->GetEntries(); i++) {
    const PFJet *pfJet = dynamic_cast<const PFJet*>(fJets->At(i));
    fPFJets->Add(pfJet);
  }

  fXlMet->Trim();
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlMet::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  ReqEventObject(fJetsName, fJets, fJetsFromBranch);
  ReqEventObject(fMuonsName, fMuons, fMuonsFromBranch);
  ReqEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  ReqEventObject(fTausName, fPFTaus, fTausFromBranch);
  ReqEventObject(fPhotonsName, fPhotons, fPhotonsFromBranch);
  ReqEventObject(fPFCandidatesName,  fPFCandidates,  fPFCandidatesFromBranch);
  ReqEventObject(fPVName, fPV, fPVFromBranch);
  ReqEventObject(fRawMetName, fRawMet, true);
  ReqEventObject(fPileUpDenName, fPileUpDen, true);

  // Create the new output collection
  fXlMet = new XlMetArr(16,fXlMetName);
  // Publish collection for further usage in the analysis
  if (fPublishOutput)
    PublishObj(fXlMet);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlMet::SlaveTerminate()
{

}
