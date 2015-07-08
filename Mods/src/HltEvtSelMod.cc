#include "MitMonoJet/Mods/interface/HltEvtSelMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
// Root headers
#include <TDataMember.h>
#include <TTree.h>
#include <TFile.h>

using namespace mithep;
ClassImp(mithep::HltEvtSelMod)

//--------------------------------------------------------------------------------------------------
HltEvtSelMod::HltEvtSelMod(const char *name, const char *title) : 
  BaseMod         (name,title),
  fEvtHdrName     (Names::gkEvtHeaderBrn),
  fEventHeader    (0),
  fVertexName     (Names::gkPVBeamSpotBrn),
  fVertices       (0),
  fHltObjsName    (""),
  fMetName        ("PFMet"),
  fMet            (0),
  fCleanJetsName  (""),
  fTrgMatchDeltaR (0.5),
  fPtJetCut       (0.),
  fEtaJetCut      (0.),
  fTupleName      ("hMonoJetTreeHlt"),
  fMonoJetEventHlt(0),
  fMonoJetTupleHlt(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HltEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.  
}

//--------------------------------------------------------------------------------------------------
void HltEvtSelMod::Process()
{
  // Process entries of the tree

  // ---- Loading what we need to proceed

  // Load the Header
  LoadBranch(fEvtHdrName);

  // Load the vertices (hopefully created before)
  fVertices = GetObject<VertexOArr>(fVertexName);
  
  // Load the HLT objects (hopefully created before)
  const TriggerObjectCol *hltObjs = GetHLTObjects(fHltObjsName);
  if (!hltObjs){
    printf(" HltEvtSelMod::Process -- TriggerObjectCol not found\n");
    // no point continuing as trigger matching will fail
    return;
  }

  // Get the MET branch
  LoadEventObject(fMetName,fMet,true);

  // Obtain all the good jets from the event cleaning module
  JetOArr *cleanJets = GetObject<JetOArr>(fCleanJetsName);

  // ---- Filling ntuple variables

  // From the top
  fMonoJetEventHlt->fRunNumber   = fEventHeader->RunNum();
  fMonoJetEventHlt->fNumVertices = fVertices->GetEntries();
  fMonoJetEventHlt->fNjets       = cleanJets->GetEntries();
  fMonoJetEventHlt->fMetEt       = fMet->At(0)->Pt();

  // Prepare all the bits necessary to get the HLT-offline match information
  fMonoJetEventHlt->ResetTriggerBits();

  // Loop on trigger objects and set the bits
  for (unsigned int iHlt=0; iHlt<hltObjs->GetEntries(); ++iHlt) {

    const TriggerObject* to = hltObjs->At(iHlt);
      
    // Muon Hlt
    if (to->Type() == TriggerObject::TriggerMuon &&
	TString(to->TrigName()).Contains("IsoMu24_eta2p1_v"))
      fMonoJetEventHlt->fHltMu = true;

    // Jet HLT matching, loop on clean jets and check if they are matched to the relevant HLT jets
    for (unsigned int iJet = 0; iJet < cleanJets->GetEntries(); iJet++) {
      // Make: object is jet and it matche within deltaR < fTrgMatchDeltaR the offline jet
      if (to->Type() == TriggerObject::TriggerJet &&
	  MathUtils::DeltaR(cleanJets->At(iJet)->Mom(),to->Mom())<0.5) {
	// Check for the two relevant bits
        if (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v"))
	  fMonoJetEventHlt->fHltJets[iJet] = true;
        if (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v"))
	  fMonoJetEventHlt->fHltJets[iJet] = true;
      }
      // Also store sme masic parameters about the jet
      fMonoJetEventHlt->fJetPt [iJet] = cleanJets->At(iJet)->Mom().Pt();
      fMonoJetEventHlt->fJetEta[iJet] = cleanJets->At(iJet)->Mom().Eta();
    }
    
    // Was there any jet trigger matching our bits? no matched offline jet required
    if (to->Type() == TriggerObject::TriggerJet) {
      if (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v") ||
	  TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v")  )
	fMonoJetEventHlt->fHltJetsNoMatch = true;
    }

    // Fill the MET HLT bits
    if (to->Type() == TriggerObject::TriggerMET &&
	TString(to->TrigName()).Contains("HLT_MET120_HBHENoiseCleaned_v"))
      fMonoJetEventHlt->fHltMets[0] = true;
    if (to->Type() == TriggerObject::TriggerMHT && 
        (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v") ||
	 TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v")) )
      fMonoJetEventHlt->fHltMets[1] = true;
  }

  // ---- Commit all to the container

  fMonoJetTupleHlt->Fill();
}
//--------------------------------------------------------------------------------------------------
void HltEvtSelMod::SlaveBegin()
{
  // Request all the things we need later on

  ReqEventObject(fMetName,   fMet, true);
  ReqBranch     (fEvtHdrName,fEventHeader);


  // ---- Create Ntuple Tree  

  fMonoJetEventHlt = new MonoJetEventHlt;
  TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  fMonoJetTupleHlt = new TTree(fTupleName.Data(),fTupleName.Data());
    

  // ---- Wow, this is ugly!! and for what? just because you do not like dictionaries?


  // Make flattish tree from classes so we do not have to rely on dictionaries for reading later
  TClass *eclass = TClass::GetClass("mithep::MonoJetEventHlt");
  TList  *elist  = eclass->GetListOfDataMembers();

  for (int i=0; i<elist->GetEntries(); ++i) {

    const TDataMember *tdm = static_cast<const TDataMember*>(elist->At(i)); //ming

    if (! (tdm->IsBasic() && tdm->IsPersistent()))
      continue;
    
    // Determine the datatype string to be used
    TString ts = TypeString(tdm);
    if (ts = TString(""))
      continue;

    // Find the address of our object
    Char_t *addr = (Char_t*) fMonoJetEventHlt; // ming:?
    assert(sizeof(Char_t) == 1);               // why do we need this?

    if      (TString(tdm->GetName()).CompareTo("fHltJets") == 0)
      fMonoJetTupleHlt->Branch(tdm->GetName(),addr+tdm->GetOffset(),
			       TString::Format("%s[4]/%s",tdm->GetName(),ts.Data())); 
    else if (TString(tdm->GetName()).CompareTo("fHltMets") == 0)
      fMonoJetTupleHlt->Branch(tdm->GetName(),addr+tdm->GetOffset(),
			       TString::Format("%s[2]/%s",tdm->GetName(),ts.Data())); 
    else if (TString(tdm->GetName()).CompareTo("fJetPt") == 0)
      fMonoJetTupleHlt->Branch(tdm->GetName(),addr+tdm->GetOffset(),
			       TString::Format("%s[4]/%s",tdm->GetName(),ts.Data())); 
    else if (TString(tdm->GetName()).CompareTo("fJetEta") == 0)
      fMonoJetTupleHlt->Branch(tdm->GetName(),addr+tdm->GetOffset(),
			       TString::Format("%s[4]/%s",tdm->GetName(),ts.Data())); 
    else
      fMonoJetTupleHlt->Branch(tdm->GetName(),addr+tdm->GetOffset(),
			       TString::Format("%s/%s",tdm->GetName(),ts.Data()));
  }
  
  AddOutput(fMonoJetTupleHlt);
}

//--------------------------------------------------------------------------------------------------
void HltEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void HltEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
TString HltEvtSelMod::TypeString(const TDataMember *tdm)
{
  TString ts= "";
  if      (TString(tdm->GetTypeName()).BeginsWith("Char_t"))    ts = "B";
  else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t"))   ts = "b";
  else if (TString(tdm->GetTypeName()).BeginsWith("Short_t"))   ts = "S";
  else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t"))  ts = "s";
  else if (TString(tdm->GetTypeName()).BeginsWith("Int_t"))     ts = "I";
  else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t"))    ts = "i";
  else if (TString(tdm->GetTypeName()).BeginsWith("Float_t"))   ts = "F";
  else if (TString(tdm->GetTypeName()).BeginsWith("Double_t"))  ts = "D";
  else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t"))  ts = "L";
  else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) ts = "l";
  else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t"))    ts = "O";
  return ts;
}
