#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/TreeFiller/interface/FillerXsIsoParticles.h"

#include "MitAna/DataTree/interface/XsIsoParticle.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;

ClassImp(mithep::FillerXsIsoParticles)

//--------------------------------------------------------------------------------------------------
FillerXsIsoParticles::FillerXsIsoParticles(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fFillXsMuons (kTRUE),
  fFillXsElectrons (kTRUE),
  fFillXsTaus (kTRUE),
  fFillXsPhotons (kTRUE),
  fPublishOutput (kTRUE),
  fMuonsName (Names::gkMuonBrn),
  fMuonsFromBranch (kTRUE),
  fMuons (0),
  fIsoMuonsName (Names::gkMuonBrn),
  fIsoMuonsFromBranch (kTRUE),
  fIsoMuons (0),
  fElectronsName (Names::gkElectronBrn),
  fElectronsFromBranch (kTRUE),
  fElectrons (0),
  fTausName (Names::gkPFTauBrn),
  fTausFromBranch (kTRUE),
  fTaus (0),
  fPhotonsName (Names::gkPhotonBrn),
  fPhotonsFromBranch (kTRUE),
  fPhotons (0),  
  fXsMuonsName ("XsMuons"),
  fXsElectronsName ("XsElectrons"),
  fXsTausName ("XsTaus"),
  fXsPhotonsName ("XsPhotons")
{
  // Constructor.
}

FillerXsIsoParticles::~FillerXsIsoParticles()
{
  // Destructor
  if (fXsMuons)
    delete fXsMuons;
  if (fXsElectrons)
    delete fXsElectrons;
  if (fXsTaus)
    delete fXsTaus;
  if (fXsPhotons)
    delete fXsPhotons;
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::Process()
{
  // make sure the out collections are empty before starting
  fXsMuons->Delete();  
  fXsElectrons->Delete();  
  fXsTaus->Delete();  
  fXsPhotons->Delete();  
  
  // Load the branches we want to work with

  if (fFillXsMuons) {
    fMuons = GetObject<MuonOArr>(fMuonsName);
    fIsoMuons = GetObject<MuonOArr>(fIsoMuonsName);
    // LoadEventObject(fMuonsName,fMuons,fMuonsFromBranch);
    // LoadEventObject(fIsoMuonsName,fIsoMuons,fIsoMuonsFromBranch);
  }
  if (fFillXsElectrons) 
    fElectrons = GetObject<ElectronOArr>(fElectronsName);
    // LoadEventObject(fElectronsName,fElectrons,fElectronsFromBranch);
  if (fFillXsTaus) 
    fTaus = GetObject<PFTauOArr>(fTausName);
    // LoadEventObject(fTausName,fTaus,fTausFromBranch);
  if (fFillXsPhotons) 
    fPhotons = GetObject<PhotonOArr>(fPhotonsName);
    // LoadEventObject(fPhotonsName,fPhotons,fPhotonsFromBranch);

  // Muon block
  if (fFillXsMuons)
    for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
      const Particle *muon = dynamic_cast<const Particle*>(fMuons->At(i));
      // Determine if the muon is passing tight selections
      Bool_t isTight = IsTightMuon(fMuons->At(i));
      // Determine if the muon is passing isolation requirements
      Bool_t isIso = IsIsoMuon(fMuons->At(i));
      
      // Fill the XsMuons collection with the reduced muon object
      FillXsIsoParticle(fXsMuons,muon,isTight,isIso);            
   }

  // Electrons block
  if (fFillXsElectrons)
    for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
      const Particle *electron = dynamic_cast<const Particle*>(fElectrons->At(i));
      // Fill the XsElectrons collection with the reduced electron object
      FillXsIsoParticle(fXsElectrons,electron);            
    }

  // Taus block
  if (fFillXsTaus)
    for (UInt_t i=0; i<fTaus->GetEntries(); ++i) {
      const Particle *tau = dynamic_cast<const Particle*>(fTaus->At(i));
      // Fill the XsTaus collection with the reduced electron object
      FillXsIsoParticle(fXsTaus,tau);            
    }

  // Photons block
  if (fFillXsPhotons)
    for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
      const Particle *photon = dynamic_cast<const Particle*>(fPhotons->At(i));
      // Fill the XsPhotons collection with the reduced electron object
      FillXsIsoParticle(fXsPhotons,photon);            
    }

  // Trim the output collections
  if (fFillXsMuons)
    fXsMuons->Trim();  
  if (fFillXsElectrons) 
    fXsElectrons->Trim();  
  if (fFillXsTaus) 
    fXsTaus->Trim();  
  if (fFillXsPhotons) 
    fXsPhotons->Trim();  
    
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. 
  if (fFillXsMuons) {
    ReqEventObject(fMuonsName,fMuons,fMuonsFromBranch);
    ReqEventObject(fIsoMuonsName,fIsoMuons,fIsoMuonsFromBranch);
  }
  if (fFillXsElectrons) 
    ReqEventObject(fElectronsName,fElectrons,fElectronsFromBranch);
  if (fFillXsTaus) 
    ReqEventObject(fTausName,fTaus,fTausFromBranch);
  if (fFillXsPhotons) 
    ReqEventObject(fPhotonsName,fPhotons,fPhotonsFromBranch);

  // Create the new output collection
  fXsMuons = new XsIsoParticleArr(16,fXsMuonsName);
  fXsElectrons = new XsIsoParticleArr(16,fXsElectronsName);
  fXsTaus = new XsIsoParticleArr(16,fXsTausName);
  fXsPhotons = new XsIsoParticleArr(16,fXsPhotonsName);
  // Publish collection for further usage in the analysis
  if (fPublishOutput) {
    PublishObj(fXsMuons);
    PublishObj(fXsElectrons);
    PublishObj(fXsTaus);
    PublishObj(fXsPhotons);
  }
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::FillXsIsoParticle(XsIsoParticleArr *pXsArr, const Particle *pParticle)
{
  // Prepare and store in an array a new XsIsoParticle 
  XsIsoParticle *thisXsIsoParticle = pXsArr->Allocate();
  new (thisXsIsoParticle) XsIsoParticle();

  // Determine basic particle proprieties
  thisXsIsoParticle->SetMom(pParticle->Mom());
  thisXsIsoParticle->SetCharge(pParticle->Charge());
    
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::FillXsIsoParticle(XsIsoParticleArr *pXsArr, const Particle *pParticle,
                                             Bool_t isTight, Bool_t isIso)
{
  // Prepare and store in an array a new XsIsoParticle 
  XsIsoParticle *thisXsIsoParticle = pXsArr->Allocate();
  new (thisXsIsoParticle) XsIsoParticle();
  
  // Determine basic particle proprieties
  thisXsIsoParticle->SetMom(pParticle->Mom());
  thisXsIsoParticle->SetCharge(pParticle->Charge());
    
  // Determine particle quality
  if (isTight) {
    if (isIso)
      thisXsIsoParticle->SetParticleId(XsIsoParticle::EParticleId::eIsoMuon);
    else 
      thisXsIsoParticle->SetParticleId(XsIsoParticle::EParticleId::eTightMuon);
  }
       
  return;
}

//--------------------------------------------------------------------------------------------------
Bool_t FillerXsIsoParticles::IsTightMuon(const Muon *muon)
{
  Bool_t theDecision = false;

  theDecision = 
  ((muon->HasGlobalTrk() && muon->GlobalTrk()->Chi2()/muon->GlobalTrk()->Ndof() < 10 
  && (muon->NSegments() > 1 || muon->NMatches() > 1) && muon->NValidHits() > 0 ) 
  || muon->IsTrackerMuon() ) &&
  (muon->BestTrk() != 0 && muon->BestTrk()->NHits() > 10 && 
  (muon->NSegments() > 1 || muon->NMatches() > 1) && muon->BestTrk()->NPixelHits() > 0 );
  
  return theDecision;
}

//--------------------------------------------------------------------------------------------------
Bool_t FillerXsIsoParticles::IsIsoMuon(const Muon *muon)
{
  float minDr = 999.;
  float deltaR = 0.01;

  // Loop on the isolated muon collection 
  // and find a matched object. Matching radius is 0.01
  for (UInt_t i=0; i<fIsoMuons->GetEntries(); ++i) {
    const Muon *isomuon = fIsoMuons->At(i);
    // compute the dR
    float thisDr = MathUtils::DeltaR(*muon, *isomuon);
    if (thisDr < minDr) 
      minDr = thisDr;
    // check if user match condition is met
    if (minDr < deltaR)
      return true;
  } // end loop on trigger objects

  return false;
}
