#ifndef MitSynchTree_H
#define MitSynchTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

class MitSynchTree {
 public:

  /// variables
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  float          met_;
  Float_t        phopt_;
  Float_t        phoeta_;
  Float_t        jetpt_;
  Float_t        jeteta_;
  
 public:
  /// this is the main element
  TFile *f_;
  TTree *tree_;

  /// default constructor  
  MitSynchTree() {}
  /// default destructor
  ~MitSynchTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitSynchTree
  void LoadTree(const char* file, const char* treeName){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(treeName));
    assert(tree_);
  }

  /// create a MitSynchTree
  void CreateTree(const char* treeName){
    tree_ = new TTree(treeName,"Limit tree");
    f_ = 0;
    InitVariables();
    
    //book the branches    
    tree_->Branch("event"           , &event_           ,   "event/i");
    tree_->Branch("run"             , &run_             ,   "run/i");
    tree_->Branch("lumi"            , &lumi_            ,   "lumi/i");
    tree_->Branch("met_"            , &met_             ,   "met_/F");
    tree_->Branch("phopt_"          , &phopt_           ,   "phopt_/F");
    tree_->Branch("phoeta_"         , &phoeta_          ,   "phoeta_/F");
    tree_->Branch("jetpt_"          , &jetpt_           ,   "jetpt_/F");
    tree_->Branch("jeteta_"         , &jeteta_          ,   "jeteta_/F");
  }

  /// create a MitSynchTree
  void DestroyTree(){
    tree_->Delete();
  }

  // initialze a MitSynchTree
  void InitTree(){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("event"           , &event_  );
    tree_->SetBranchAddress("run"             , &run_    );
    tree_->SetBranchAddress("lumi"            , &lumi_   );
    tree_->SetBranchAddress("met_"            , &met_    );
    tree_->SetBranchAddress("phopt_"          , &phopt_  );
    tree_->SetBranchAddress("phoeta_"         , &phoeta_ );
    tree_->SetBranchAddress("jetpt_"          , &jetpt_  );
    tree_->SetBranchAddress("jeteta_"         , &jeteta_ );

    gErrorIgnoreLevel = currentState;
  }

}; 

inline void 
MitSynchTree::InitVariables(){
  // inizialize variables
  event_         = 0;
  run_           = 0;
  lumi_          = 0;
  met_           = -.1;
  phopt_         = -.1;
  phoeta_        = -.1;
  jetpt_         = -.1;
  jeteta_        = -.1;

}

#endif
