#ifndef MitLimitTree_H
#define MitLimitTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

class MitLimitTree {
 public:

  /// variables
  unsigned int     event_;
  unsigned int     run_;
  unsigned int     lumi_;
  Float_t          mvamet_;
  Float_t          mvametphi_;
  Float_t          mt_;
  Float_t          mll_;
  Float_t          ptll_;
  Float_t          ptpho_;
  Int_t            njets_;
  Float_t          jet1pt_;
  Float_t          genjetpt_;
  Float_t          genVpt_;
  Float_t          genVphi_;
  Float_t          dmpt_;
  Float_t          weight_;
  
 public:
  /// this is the main element
  TFile *f_;
  TTree *tree_;

  /// default constructor  
  MitLimitTree() {}
  /// default destructor
  ~MitLimitTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitLimitTree
  void LoadTree(const char* file, const char* treeName){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(treeName));
    assert(tree_);
  }

  /// create a MitLimitTree
  void CreateTree(const char* treeName){
    tree_ = new TTree(treeName,"Limit tree");
    f_ = 0;
    InitVariables();
    
    //book the branches    
    tree_->Branch("event"          , &event_           ,   "event/i");
    tree_->Branch("run"            , &run_             ,   "run/i");
    tree_->Branch("lumi"           , &lumi_            ,   "lumi/i");
    tree_->Branch("mvamet"         , &mvamet_          ,   "mvamet/F");
    tree_->Branch("mvametphi"      , &mvametphi_       ,   "mvametphi/F");
    tree_->Branch("mt"             , &mt_              ,   "mt/F");
    tree_->Branch("mll"            , &mll_             ,   "mll/F");
    tree_->Branch("ptll"           , &ptll_            ,   "ptll/F");
    tree_->Branch("ptpho"          , &ptpho_           ,   "ptpho/F");
    tree_->Branch("njets"          , &njets_           ,   "njets/i");
    tree_->Branch("jet1pt"         , &jet1pt_          ,   "jet1pt/F");
    tree_->Branch("genjetpt"       , &genjetpt_        ,   "genjetpt/F");
    tree_->Branch("genVpt"         , &genVpt_          ,   "genVpt/F");
    tree_->Branch("genVphi"        , &genVphi_         ,   "genVphi/F");
    tree_->Branch("dmpt"           , &dmpt_            ,   "dmpt/F");
    tree_->Branch("weight"         , &weight_          ,   "weight/F");
  }

  /// create a MitLimitTree
  void DestroyTree(){
    tree_->Delete();
  }

  // initialze a MitLimitTree
  void InitTree(){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("event"          , &event_);
    tree_->SetBranchAddress("run"            , &run_);
    tree_->SetBranchAddress("lumi"           , &lumi_);
    tree_->SetBranchAddress("mvamet"         , &mvamet_);
    tree_->SetBranchAddress("mvametphi"      , &mvametphi_);
    tree_->SetBranchAddress("mt"             , &mt_);
    tree_->SetBranchAddress("mll"            , &mll_);
    tree_->SetBranchAddress("ptll"           , &ptll_);
    tree_->SetBranchAddress("ptpho"          , &ptpho_);
    tree_->SetBranchAddress("njets"          , &njets_);
    tree_->SetBranchAddress("jet1pt"         , &jet1pt_);
    tree_->SetBranchAddress("genjetpt"       , &genjetpt_);
    tree_->SetBranchAddress("genVpt"         , &genVpt_);
    tree_->SetBranchAddress("genVphi"        , &genVphi_);
    tree_->SetBranchAddress("dmpt"           , &dmpt_);
    tree_->SetBranchAddress("weight"         , &weight_);

    gErrorIgnoreLevel = currentState;
  }

}; 

inline void 
MitLimitTree::InitVariables(){
  // inizialize variables
  event_     = 0;
  run_       = 0;
  lumi_      = 0;
  mvamet_    = -1.;
  mvametphi_ = -1.;
  mt_        = -1;
  mll_       = -1;
  ptll_      = -1;
  ptpho_     = -1;
  njets_     = -1;
  jet1pt_    = -1.;
  genjetpt_  = -1.;
  genVpt_    = -1.;
  genVphi_   = -1.;
  dmpt_      = -1.;
  weight_    = -1.;

}

#endif
