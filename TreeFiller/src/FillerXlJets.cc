#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlJets.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitAna/DataTree/interface/XlJet.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;

ClassImp(mithep::FillerXlJets)

//--------------------------------------------------------------------------------------------------
FillerXlJets::FillerXlJets(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fQGTaggingActive (kTRUE),
  fQGTaggerCHS (kFALSE),
  fPublishOutput (kTRUE),
  fJetsName (Names::gkPFJetBrn),
  fJetsFromBranch (kTRUE),
  fJets (0),
  fPileUpDenName(Names::gkPileupEnergyDensityBrn),
  fPileUpDenFromBranch(kTRUE),
  fPileUpDen(0),
  fVertexesName (ModNames::gkGoodVertexesName),
  fVertexesFromBranch(kFALSE),
  fVertexes(0),
  fXlJetsName ("XlJets")
{
  // Constructor.
}

FillerXlJets::~FillerXlJets()
{
  // Destructor
  if (fXlJets)
    delete fXlJets;
  
  delete fQGTagger;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::Process()
{
  // make sure the out collections are empty before starting
  fXlJets->Delete();  
  
  // // Load the branches we want to work with
  // LoadEventObject(fJetsName,fJets,fJetsFromBranch);
  // if (fQGTaggingActive) {
  //   LoadEventObject(fPileUpDenName,fPileUpDen,fPileUpDenFromBranch);
  //   LoadEventObject(fVertexesName,fVertexes,fVertexesFromBranch);
  // }
 
  fJets = GetObject<JetOArr>(fJetsName);
  if (fQGTaggingActive){
    fPileUpDen = GetObject<PileupEnergyDensityCol>(fPileUpDenName);
    fVertexes = GetObject<VertexCol>(fVertexesName);
  }

  // Setup pileup density for QG computation
  if (fQGTaggingActive)
    fQGTagger->SetRhoIso(fPileUpDen->At(0)->RhoRandomLowEta());
  
  // Loop over jets
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {
      
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(i));
    if (! jet) {
      printf(" FillerXlJets::Process() - ERROR - jets provided are not PFJets.");
      break;
    }
 
    // mark jet (and consequently its consituents) for further use in skim
    jet->Mark();       
    
    // perform jet analysis and fill the extended XlJet object
    FillXlJet(jet);      
    
  }    
  // Trim the output collection
  fXlJets->Trim();
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. 
  ReqEventObject(fJetsName,fJets,fJetsFromBranch);
  ReqEventObject(fPileUpDenName,fPileUpDen,fPileUpDenFromBranch);
  ReqEventObject(fVertexesName,fVertexes,fVertexesFromBranch);

  // Initialize area caculation (done with ghost particles)

  // Create the new output collection
  fXlJets = new XlJetArr(16,fXlJetsName);
  // Publish collection for further usage in the analysis
  if (fPublishOutput)
    PublishObj(fXlJets);
  
  // Initialize QGTagger class
  fQGTagger = new QGTagger(fQGTaggerCHS);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlJet(const PFJet *pPFJet)
{
  // Prepare and store in an array a new FatJet 
  XlJet *newXlJet = fXlJets->Allocate();
  new (newXlJet) XlJet(*pPFJet);

  // Compute and store weighted charge
  newXlJet->SetCharge();
  
  // Prepare and store QG tagging info
  float qgValue = -1.;
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(pPFJet, fVertexes);
    qgValue = fQGTagger->QGValue();
  }
  newXlJet->SetQGTag(qgValue);
    
  // Prepare and store jet pull info
  TVector2 newXlJetPull = GetPull(pPFJet);
  newXlJet->SetPullY(newXlJetPull.X());
  newXlJet->SetPullPhi(newXlJetPull.Y());
      
  return;
}

//--------------------------------------------------------------------------------------------------
TVector2 FillerXlJets::GetPull(const PFJet *inPFJet)
{
  double dYSum   = 0;
  double dPhiSum = 0;
  const unsigned int nPFCands = inPFJet->NPFCands();

  // Loop on input jet constituents and get the color pull  
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {    
    const PFCandidate *pfCand = inPFJet->PFCand(ipf);
    double pt_i=0, y_i=0, phi_i=0;
    pt_i = pfCand->Pt();
    y_i = pfCand->Rapidity();
    phi_i = pfCand->Phi();

    double dY   = y_i - inPFJet->Rapidity();
    double dPhi = MathUtils::DeltaPhi(phi_i,inPFJet->Phi());
    double weight = pt_i*sqrt(dY*dY + dPhi*dPhi);
    dYSum += weight*dY;
    dPhiSum += weight*dPhi;
  }

  return TVector2(dYSum/inPFJet->Pt(), dPhiSum/inPFJet->Pt());    
}
