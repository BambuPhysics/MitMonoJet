// $Id: $

#include "MitMonoJet/Mods/interface/SkimJetsMod.h"

using namespace mithep;

ClassImp(mithep::SkimJetsMod)

//--------------------------------------------------------------------------------------------------
SkimJetsMod::SkimJetsMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fBranchName("SkmSetMe"),
  fColFromBranch(kTRUE),
  fColMarkFilter(kTRUE),
  fPublishArray(kFALSE),
  fCol(0),
  fColSkm(0),
  fArrSkm(0)
{
  // Constructor.
}

SkimJetsMod::~SkimJetsMod() 
{
  // Destructor.
  if (fColSkm)
    delete fColSkm;

  if (fArrSkm)
    delete fArrSkm;
}

//--------------------------------------------------------------------------------------------------
void SkimJetsMod::Process()
{
  // make sure the collection is empty before starting
  if (!fPublishArray)
    fColSkm->Reset();
  else
    fArrSkm->Reset();

  // load the branch with the proper method
  if (fColFromBranch)
    LoadBranch(GetBranchName());
  else 
    fCol = GetObject<JetCol>(GetBranchName());
    // LoadEventObject(GetBranchName(), fCol);

  // loop on the input collection and apply the filter on mark if required
  const UInt_t entries = fCol->GetEntries();
  for (UInt_t i=0; i<entries; ++i) {

    // if marked or all objects requested
    if (!fColMarkFilter || fCol->At(i)->IsMarked()) {

      // convert the input jet into a pfjet
      const PFJet* inpfjet = dynamic_cast<const PFJet*>(fCol->At(i));
      // Make sure the mark is removed to avoid having the future objects already marked
      // - collections that are not skimmed will have dangeling marks but that us fine
      //   as marking is exclusively used by the skimmer
      if (fColMarkFilter)
        fCol->At(i)->UnmarkMe();

      // fill the output Array/ObjArray
      if (!fPublishArray)
        fColSkm->Add(inpfjet);
      else {        
        PFJet *outpfjet = fArrSkm->Allocate();
        new (outpfjet) PFJet(*inpfjet);
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
void SkimJetsMod::SlaveBegin()
{
  // Request the marked input branch
  if (fColFromBranch)
    ReqBranch(GetBranchName(), fCol);
  else 
    ReqEventObject(GetBranchName(), fCol, fColFromBranch);

  // Request the branch to be published
  if (!fPublishArray) {
    fColSkm = new ObjArray<PFJet>(0,TString("Skm")+GetBranchName());
    PublishObj(fColSkm);
  }
  else {
    fArrSkm = new Array<PFJet>(0,TString("Skm")+GetBranchName());    
    PublishObj(fArrSkm);
  }
}

//--------------------------------------------------------------------------------------------------
void SkimJetsMod::SlaveTerminate()
{
}
