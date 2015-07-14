#ifndef MitGPTree_H
#define MitGPTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

class MitGPTree {
 public:
  /// float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  /// DON'T CHANGE ORDER
  enum DataType {
    data,
    dmg,
    dmj,
    top,
    dyll,
    wjets,
    vv,
    wgamma,
    vvv,
    gj,
    qcd,
    other
  };

  enum Selection {
    GoodPhoton      = 1UL<<0,    // good photon
    Fake            = 1UL<<1,    // fake
    BeamHalo        = 1UL<<2,    // beam halo
    DiLepton        = 1UL<<3    // dilepton
  };

  /// variables
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   trigger_;
  unsigned int   HLTmatch_;
  unsigned int   nvtx_;
  unsigned int   cuts_;
  float          scale1fb_;
  float          genmet_;
  float          genmetPhi_;
  float          metRaw_;
  float          metRawPhi_;
  float          met_;
  float          metPhi_;
  float          mvamet_;
  float          mvametPhi_;
  float          mvaCov00_;
  float          mvaCov10_;
  float          mvaCov01_;
  float          mvaCov11_;
  float          metCorZ_;
  float          metCorZPhi_;
  float          metCorW_;
  float          metCorWPhi_;
  float          metRawCorZ_;
  float          metRawCorZPhi_;
  float          metRawCorW_;
  float          metRawCorWPhi_;
  float          mvametCorZ_;
  float          mvametCorZPhi_;
  float          mvametCorW_;
  float          mvametCorWPhi_;
  float          sumEt_;
  float          metSig_;
  DataType       dstype_;

  LorentzVector  genZ_;
  LorentzVector  genH_;
  LorentzVector  genMuon1_;
  LorentzVector  genMuon2_;

  unsigned int   nlep_;
  LorentzVector  lep1_;
  int            lid1_;
  int		 lep1IsTightMuon_;
  int		 lep1IsIsolated_;
  float          lep1PtErr_;
  LorentzVector  lep2_;
  int            lid2_;
  int		 lep2IsTightMuon_;
  int		 lep2IsIsolated_;
  float          lep2PtErr_;
  LorentzVector  lep3_;
  int            lid3_;
  int		 lep3IsTightMuon_;

  unsigned int   ntaus_;
  LorentzVector  tau1_;
  LorentzVector  tau2_;

  unsigned int   nphotons_;
  LorentzVector  pho1_;
  LorentzVector  pho2_;
/*   LorentzVector  pho3_; */
/*   LorentzVector  pho4_; */

  unsigned int   njets_;
  unsigned int   noiseCleaning_;
  LorentzVector  jet1_;
  float          jet1CHF_;
  float          jet1NHF_;
  float          jet1NEMF_;
  float          jet1Unc_;
  float          jet1Btag_;
  float          jet1QGtag_;
  unsigned int   jet1PartonId_;
  float          jet1QGRho_;
  float          jet1QGPtD_;
  float          jet1QGAxis1_;
  float          jet1QGAxis2_;
  float          jet1QGMult_;
  LorentzVector  jet2_;
  float          jet2CHF_;
  float          jet2NHF_;
  float          jet2NEMF_;
  float          jet2Unc_;
  float          jet2Btag_;
  float          jet2QGtag_;
  LorentzVector  jet3_;
  float          jet3CHF_;
  float          jet3NHF_;
  float          jet3NEMF_;
  float          jet3Btag_;
  LorentzVector  jet4_;
  float          jet4Btag_;

  unsigned int   ntracks_;
  LorentzVector  track1_;
  LorentzVector  track2_;
  LorentzVector  track3_;

  float          Q_;
  float          id1_;
  float          x1_;
  float          pdf1_;
  float          id2_;
  float          x2_;
  float          pdf2_;
  int            processId_;
  float          npu_;
  float          npuPlusOne_;
  float          npuMinusOne_;
  int            metFiltersWord_;

 public:
  /// this is the main element
  TFile *f_;
  TTree *tree_;

  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  MitGPTree():
    genZPtr_(&genZ_),
    genHPtr_(&genH_),
    genMuon1Ptr_(&genMuon1_),genMuon2Ptr_(&genMuon2_),
    lepPtr1_(&lep1_),lepPtr2_(&lep2_),lepPtr3_(&lep3_),
    tauPtr1_(&tau1_),tauPtr2_(&tau2_),
    phoPtr1_(&pho1_),phoPtr2_(&pho2_),//phoPtr3_(&pho3_),phoPtr4_(&pho4_),
    jetPtr1_(&jet1_),jetPtr2_(&jet2_),jetPtr3_(&jet3_),jetPtr4_(&jet4_),
    trackPtr1_(&track1_),trackPtr2_(&track2_),trackPtr3_(&track3_){}
  /// default destructor
  ~MitGPTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitGPTree
  void LoadTree(const char* file, int type = -1){
    // to load three different ntuples in the same job DMTree0/1
    // type == 0/1 if all variables was added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    f_ = TFile::Open(file);
    assert(f_);
    if     (type == 0) tree_ = dynamic_cast<TTree*>(f_->Get("MJetTree"));
    else if(type == 1) tree_ = dynamic_cast<TTree*>(f_->Get("MJetTreeDiLepton"));
    else if(type == 2) tree_ = dynamic_cast<TTree*>(f_->Get("MJetTreeWlnu"));
    else	       tree_ = dynamic_cast<TTree*>(f_->Get("tree"));
    assert(tree_);
  }

  /// load a MitGPTree
  void LoadTree(const char* file, const char* treeName){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get(treeName));
    assert(tree_);
  }

  /// create a MitGPTree
  void CreateTree(int type = -1){
    assert(type==type); // just to suppress warnings
    // to create three different ntuples in the same job DMTree0/1
    // type == 0/1 add all variables
    // type = -1 (default) add a minimum set of variables with tree as name
    if     (type == 0) tree_ = new TTree("MJetTree","Smurf ntuples");
    else if(type == 1) tree_ = new TTree("MJetTreeDiLepton","Smurf ntuples");
    else if(type == 2) tree_ = new TTree("MJetTreeWlnu","Smurf ntuples");
    else               tree_ = new TTree("tree","Smurf ntuples");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("event"        , &event_        ,   "event/i");
    tree_->Branch("run"          , &run_          ,   "run/i");
    tree_->Branch("lumi"         , &lumi_         ,   "lumi/i");
    tree_->Branch("trigger"      , &trigger_      ,   "trigger/i");
    tree_->Branch("HLTmatch"     , &HLTmatch_     ,   "HLTmatch/i");
    tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/i");
    tree_->Branch("scale1fb"     , &scale1fb_     ,   "scale1fb/F");
    tree_->Branch("genmet"       , &genmet_       ,   "genmet/F");
    tree_->Branch("genmetPhi"    , &genmetPhi_    ,   "genmetPhi/F");
    tree_->Branch("metRaw"       , &metRaw_       ,   "metRaw/F");
    tree_->Branch("metRawPhi"    , &metRawPhi_    ,   "metRawPhi/F");
    tree_->Branch("met"          , &met_          ,   "met/F");
    tree_->Branch("metPhi"       , &metPhi_       ,   "metPhi/F");
    tree_->Branch("mvamet"       , &mvamet_       ,   "mvamet/F");
    tree_->Branch("mvametPhi"    , &mvametPhi_    ,   "mvametPhi/F");
    tree_->Branch("mvaCov00"     , &mvaCov00_     ,   "mvaCov00/F");
    tree_->Branch("mvaCov10"     , &mvaCov10_     ,   "mvaCov10/F");
    tree_->Branch("mvaCov01"     , &mvaCov01_     ,   "mvaCov01/F");
    tree_->Branch("mvaCov11"     , &mvaCov11_     ,   "mvaCov11/F");
    tree_->Branch("metCorZ"      , &metCorZ_      ,   "metCorZ/F");
    tree_->Branch("metCorZPhi"   , &metCorZPhi_   ,   "metCorZPhi/F");
    tree_->Branch("metCorW"      , &metCorW_      ,   "metCorW/F");
    tree_->Branch("metCorWPhi"   , &metCorWPhi_   ,   "metCorWPhi/F");
    tree_->Branch("metRawCorZ"      , &metRawCorZ_      ,   "metRawCorZ/F");
    tree_->Branch("metRawCorZPhi"   , &metRawCorZPhi_   ,   "metRawCorZPhi/F");
    tree_->Branch("metRawCorW"      , &metRawCorW_      ,   "metRawCorW/F");
    tree_->Branch("metRawCorWPhi"   , &metRawCorWPhi_   ,   "metRawCorWPhi/F");
    tree_->Branch("mvametCorZ"      , &mvametCorZ_      ,   "mvametCorZ/F");
    tree_->Branch("mvametCorZPhi"   , &mvametCorZPhi_   ,   "mvametCorZPhi/F");
    tree_->Branch("mvametCorW"      , &mvametCorW_      ,   "mvametCorW/F");
    tree_->Branch("mvametCorWPhi"   , &mvametCorWPhi_   ,   "mvametCorWPhi/F");
    tree_->Branch("sumEt"        , &sumEt_        ,   "sumEt/F");
    tree_->Branch("metSig"       , &metSig_       ,   "metSig/F");
    tree_->Branch("dstype"       , &dstype_       ,   "dstype/I");

    tree_->Branch("genZ"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genZPtr_);
    tree_->Branch("genH"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genHPtr_);
    tree_->Branch("genMuon1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genMuon1Ptr_);
    tree_->Branch("genMuon2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genMuon2Ptr_);

    tree_->Branch("nlep"         , &nlep_         ,   "nlep/i");
    tree_->Branch("lep1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lid1"         , &lid1_         ,   "lid1/I");
    tree_->Branch("lep1IsTightMuon"         , &lep1IsTightMuon_         ,   "lep1IsTightMuon/I");
    tree_->Branch("lep1IsIsolated"         , &lep1IsIsolated_         ,   "lep1IsIsolated/I");
    tree_->Branch("lep1PtErr"         , &lep1PtErr_         ,   "lep1PtErr/F");
    tree_->Branch("lep2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lid2"         , &lid2_         ,   "lid2/I");
    tree_->Branch("lep2IsTightMuon"         , &lep2IsTightMuon_         ,   "lep2IsTightMuon/I");
    tree_->Branch("lep2IsIsolated"         , &lep2IsIsolated_         ,   "lep2IsIsolated/I");
    tree_->Branch("lep2PtErr"         , &lep2PtErr_         ,   "lep2PtErr/F");
    tree_->Branch("lep3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr3_);
    tree_->Branch("lid3"         , &lid3_         ,   "lid3/I");
    tree_->Branch("lep3IsTightMuon"         , &lep3IsTightMuon_         ,   "lep3IsTightMuon/I");

    tree_->Branch("ntaus"        , &ntaus_     ,   "ntaus/i");
    tree_->Branch("tau1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tauPtr1_);
    tree_->Branch("tau2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tauPtr2_);

    tree_->Branch("nphotons"     , &nphotons_     ,   "nphotons/i");
    tree_->Branch("pho1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr1_);
    tree_->Branch("pho2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr2_);
/*     tree_->Branch("pho3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr3_); */
/*     tree_->Branch("pho4"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr4_); */

    tree_->Branch("njets"        , &njets_        ,   "njets/i");
    tree_->Branch("noiseCleaning", &noiseCleaning_ ,   "noiseCleaning/i");
    tree_->Branch("jet1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet1CHF"      , &jet1CHF_      ,   "jet1CHF/F");
    tree_->Branch("jet1NHF"      , &jet1NHF_      ,   "jet1NHF/F");
    tree_->Branch("jet1NEMF"     , &jet1NEMF_      ,   "jet1NEMF/F");
    tree_->Branch("jet1Unc"      , &jet1Unc_      ,   "jet1Unc/F");
    tree_->Branch("jet1Btag"     , &jet1Btag_     ,   "jet1Btag/F");
    tree_->Branch("jet1QGtag"     , &jet1QGtag_     ,   "jet1QGtag/F");
    tree_->Branch("jet1PartonId"     , &jet1PartonId_     ,   "jet1PartonId/i");
    tree_->Branch("jet1QGRho"     , &jet1QGRho_    ,   "jet1QGRho/F");
    tree_->Branch("jet1QGPtD"     , &jet1QGPtD_    ,   "jet1QGPtD/F");
    tree_->Branch("jet1QGAxis1"   , &jet1QGAxis1_  ,   "jet1QGAxis1/F");
    tree_->Branch("jet1QGAxis2"   , &jet1QGAxis2_  ,   "jet1QGAxis2/F");
    tree_->Branch("jet1QGMult"    , &jet1QGMult_   ,   "jet1QGMult/F");
    
    tree_->Branch("jet2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet2CHF"      , &jet2CHF_      ,   "jet2CHF/F");
    tree_->Branch("jet2NHF"      , &jet2NHF_      ,   "jet2NHF/F");
    tree_->Branch("jet2NEMF"     , &jet2NEMF_      ,   "jet2NEMF/F");
    tree_->Branch("jet2Unc"      , &jet2Unc_     ,   "jet2Unc/F");
    tree_->Branch("jet2Btag"     , &jet2Btag_     ,   "jet2Btag/F");
    tree_->Branch("jet2QGtag"     , &jet2QGtag_     ,   "jet2QGtag/F");
   
    tree_->Branch("jet3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
    tree_->Branch("jet3CHF"      , &jet3CHF_      ,   "jet3CHF/F");
    tree_->Branch("jet3NHF"      , &jet3NHF_      ,   "jet3NHF/F");
    tree_->Branch("jet3NEMF"     , &jet3NEMF_      ,   "jet3NEMF/F");
    tree_->Branch("jet3Btag"     , &jet3Btag_     ,   "jet3Btag/F");
    tree_->Branch("jet4"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr4_);
    tree_->Branch("jet4Btag"     , &jet4Btag_     ,   "jet4Btag/F");

    tree_->Branch("ntracks"        , &ntracks_        ,   "ntracks/i");
    tree_->Branch("track1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &trackPtr1_);
    tree_->Branch("track2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &trackPtr2_);
    tree_->Branch("track3"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &trackPtr3_);

    tree_->Branch("Q",             &Q_	  ,     "Q/F");
    tree_->Branch("id1",           &id1_  ,     "id1/F");
    tree_->Branch("x1",            &x1_	  ,     "x1/F");
    tree_->Branch("pdf1",          &pdf1_ ,     "pdf1/F");
    tree_->Branch("id2",           &id2_  ,     "id2/F");
    tree_->Branch("x2",            &x2_	  ,     "x2/F");
    tree_->Branch("pdf2",          &pdf2_ ,     "pdf2/F");
    tree_->Branch("processId",     &processId_  , "processId/I");
    tree_->Branch("npu",           &npu_        , "npu/F");
    tree_->Branch("npuPlusOne",    &npuPlusOne_ , "npuPlusOne/F");
    tree_->Branch("npuMinusOne",   &npuMinusOne_, "npuMinusOne/F");
    tree_->Branch("metFiltersWord",   &metFiltersWord_, "metFiltersWord/I");

  }

  // initialze a MitGPTree
  void InitTree(int type = -1){
    assert(type==type); // just to suppress warnings
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("event",         &event_);
    tree_->SetBranchAddress("run",           &run_);
    tree_->SetBranchAddress("lumi",          &lumi_);
    tree_->SetBranchAddress("trigger",       &trigger_);
    tree_->SetBranchAddress("HLTmatch",      &HLTmatch_);
    tree_->SetBranchAddress("nvtx",          &nvtx_);
    tree_->SetBranchAddress("scale1fb",      &scale1fb_);
    tree_->SetBranchAddress("genmet",        &genmet_);
    tree_->SetBranchAddress("genmetPhi",     &genmetPhi_);
    tree_->SetBranchAddress("metRaw",        &metRaw_);
    tree_->SetBranchAddress("metRawPhi",     &metRawPhi_);
    tree_->SetBranchAddress("met",           &met_);
    tree_->SetBranchAddress("metPhi",        &metPhi_);
    tree_->SetBranchAddress("mvamet",        &mvamet_);
    tree_->SetBranchAddress("mvametPhi",     &mvametPhi_);
    tree_->SetBranchAddress("mvaCov00",      &mvaCov00_);
    tree_->SetBranchAddress("mvaCov10",      &mvaCov10_);
    tree_->SetBranchAddress("mvaCov01",      &mvaCov01_);
    tree_->SetBranchAddress("mvaCov11",      &mvaCov11_);
    tree_->SetBranchAddress("metCorZ",       &metCorZ_);
    tree_->SetBranchAddress("metCorZPhi",    &metCorZPhi_);
    tree_->SetBranchAddress("metCorW",       &metCorW_);
    tree_->SetBranchAddress("metCorWPhi",    &metCorWPhi_);
    tree_->SetBranchAddress("metRawCorZ",       &metRawCorZ_);
    tree_->SetBranchAddress("metRawCorZPhi",    &metRawCorZPhi_);
    tree_->SetBranchAddress("metRawCorW",       &metRawCorW_);
    tree_->SetBranchAddress("metRawCorWPhi",    &metRawCorWPhi_);
    tree_->SetBranchAddress("mvametCorZ",       &mvametCorZ_);
    tree_->SetBranchAddress("mvametCorZPhi",    &mvametCorZPhi_);
    tree_->SetBranchAddress("mvametCorW",       &mvametCorW_);
    tree_->SetBranchAddress("mvametCorWPhi",    &mvametCorWPhi_);
    tree_->SetBranchAddress("sumEt",         &sumEt_);
    tree_->SetBranchAddress("metSig",        &metSig_);
    tree_->SetBranchAddress("dstype",        &dstype_);

    tree_->SetBranchAddress("genZ",          &genZPtr_);
    tree_->SetBranchAddress("genH",          &genHPtr_);
    tree_->SetBranchAddress("genMuon1",      &genMuon1Ptr_);
    tree_->SetBranchAddress("genMuon2",      &genMuon2Ptr_);

    tree_->SetBranchAddress("nlep",          &nlep_);
    tree_->SetBranchAddress("lep1",          &lepPtr1_);
    tree_->SetBranchAddress("lid1",          &lid1_);
    tree_->SetBranchAddress("lep1IsTightMuon",&lep1IsTightMuon_);
    tree_->SetBranchAddress("lep1IsIsolated",&lep1IsIsolated_);
    tree_->SetBranchAddress("lep1PtErr",     &lep1PtErr_);
    tree_->SetBranchAddress("lep2",          &lepPtr2_);
    tree_->SetBranchAddress("lid2",          &lid2_);
    tree_->SetBranchAddress("lep2IsTightMuon",&lep2IsTightMuon_);
    tree_->SetBranchAddress("lep2IsIsolated",&lep2IsIsolated_);
    tree_->SetBranchAddress("lep2PtErr",     &lep2PtErr_);
    tree_->SetBranchAddress("lep3",          &lepPtr3_);
    tree_->SetBranchAddress("lid3",          &lid3_);
    tree_->SetBranchAddress("lep3IsTightMuon",&lep3IsTightMuon_);

    tree_->SetBranchAddress("ntaus",         &ntaus_);
    tree_->SetBranchAddress("tau1",          &tauPtr1_);
    tree_->SetBranchAddress("tau2",          &tauPtr2_);

    tree_->SetBranchAddress("nphotons"                  , &nphotons_);
    tree_->SetBranchAddress("pho1"                      , &phoPtr1_);
    tree_->SetBranchAddress("pho2"                      , &phoPtr2_);
/*     tree_->SetBranchAddress("pho3"                      , &phoPtr3_); */
/*     tree_->SetBranchAddress("pho4"                      , &phoPtr4_); */

    tree_->SetBranchAddress("njets",         &njets_);
    tree_->SetBranchAddress("noiseCleaning", &noiseCleaning_);
    tree_->SetBranchAddress("jet1",          &jetPtr1_);
    tree_->SetBranchAddress("jet1CHF",       &jet1CHF_);
    tree_->SetBranchAddress("jet1NHF",       &jet1NHF_);
    tree_->SetBranchAddress("jet1NEMF",      &jet1NEMF_);
    tree_->SetBranchAddress("jet1Unc",       &jet1Unc_);
    tree_->SetBranchAddress("jet1Btag",      &jet1Btag_);
    tree_->SetBranchAddress("jet1QGtag",     &jet1QGtag_);
    tree_->SetBranchAddress("jet1PartonId",  &jet1PartonId_);
    tree_->SetBranchAddress("jet1QGRho"    , &jet1QGRho_   );
    tree_->SetBranchAddress("jet1QGPtD"    , &jet1QGPtD_   );
    tree_->SetBranchAddress("jet1QGAxis1"  , &jet1QGAxis1_ );
    tree_->SetBranchAddress("jet1QGAxis2"  , &jet1QGAxis2_ );
    tree_->SetBranchAddress("jet1QGMult"   , &jet1QGMult_  );
    
    tree_->SetBranchAddress("jet2",          &jetPtr2_);
    tree_->SetBranchAddress("jet2CHF",       &jet2CHF_);
    tree_->SetBranchAddress("jet2NHF",       &jet2NHF_);
    tree_->SetBranchAddress("jet2NEMF",      &jet2NEMF_);
    tree_->SetBranchAddress("jet2Unc",       &jet2Unc_);
    tree_->SetBranchAddress("jet2Btag",      &jet2Btag_);
    tree_->SetBranchAddress("jet2QGtag",     &jet2QGtag_);
    tree_->SetBranchAddress("jet3",          &jetPtr3_);
    tree_->SetBranchAddress("jet3CHF",       &jet3CHF_);
    tree_->SetBranchAddress("jet3NHF",       &jet3NHF_);
    tree_->SetBranchAddress("jet3NEMF",      &jet3NEMF_);
    tree_->SetBranchAddress("jet3Btag",      &jet3Btag_);
    tree_->SetBranchAddress("jet4",          &jetPtr4_);
    tree_->SetBranchAddress("jet4Btag",      &jet4Btag_);

    tree_->SetBranchAddress("ntracks",       &ntracks_);
    tree_->SetBranchAddress("track1",        &trackPtr1_);
    tree_->SetBranchAddress("track2",        &trackPtr2_);
    tree_->SetBranchAddress("track3",        &trackPtr3_);

    tree_->SetBranchAddress("Q",	     &Q_);
    tree_->SetBranchAddress("id1",	     &id1_);
    tree_->SetBranchAddress("x1",	     &x1_);
    tree_->SetBranchAddress("pdf1",	     &pdf1_);
    tree_->SetBranchAddress("id2",	     &id2_);
    tree_->SetBranchAddress("x2",	     &x2_);
    tree_->SetBranchAddress("pdf2",	     &pdf2_);
    tree_->SetBranchAddress("processId",     &processId_);
    tree_->SetBranchAddress("npu",	     &npu_);
    tree_->SetBranchAddress("npuPlusOne",    &npuPlusOne_);
    tree_->SetBranchAddress("npuMinusOne",   &npuMinusOne_);
    tree_->SetBranchAddress("metFiltersWord",   &metFiltersWord_);

    gErrorIgnoreLevel = currentState;
  }

  private:

  LorentzVector* genZPtr_;
  LorentzVector* genHPtr_;
  LorentzVector* genMuon1Ptr_;
  LorentzVector* genMuon2Ptr_;
  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* lepPtr3_;
  LorentzVector* tauPtr1_;
  LorentzVector* tauPtr2_;
  LorentzVector* phoPtr1_;
  LorentzVector* phoPtr2_;
  LorentzVector* phoPtr3_;
  LorentzVector* phoPtr4_;
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* jetPtr3_;
  LorentzVector* jetPtr4_;
  LorentzVector* trackPtr1_;
  LorentzVector* trackPtr2_;
  LorentzVector* trackPtr3_;
}; 

inline void 
MitGPTree::InitVariables(){
  // inizialize variables
  event_         = 0;
  run_           = 0;
  lumi_          = 0;
  trigger_       = 0;
  HLTmatch_       = 0;
  nvtx_          = 0;
  scale1fb_      = 0;
  genmet_        = -999.;
  genmetPhi_     = -999.;
  metRaw_        = -999.;
  metRawPhi_     = -999.;
  met_           = -999.;
  metPhi_        = -999.;
  mvamet_        = -999.;
  mvametPhi_     = -999.;
  mvaCov00_      = -999.;
  mvaCov10_      = -999.; 
  mvaCov01_      = -999.;
  mvaCov11_      = -999.;
  metCorZ_       = -999.;
  metCorZPhi_    = -999.;
  metCorW_       = -999.;
  metCorWPhi_    = -999.;
  metRawCorZ_       = -999.;
  metRawCorZPhi_    = -999.;
  metRawCorW_       = -999.;
  metRawCorWPhi_    = -999.;
  mvametCorZ_       = -999.;
  mvametCorZPhi_    = -999.;
  mvametCorW_       = -999.;
  mvametCorWPhi_    = -999.;

  sumEt_         = -999.;
  metSig_        = -999.;
  dstype_        = data;

  genZ_       	 = LorentzVector();
  genH_       	 = LorentzVector();
  genMuon1_      = LorentzVector();
  genMuon2_      = LorentzVector();

  nlep_          = 0;
  lep1_       	 = LorentzVector();
  lid1_          = 0;
  lep1IsTightMuon_ = 0;
  lep1IsIsolated_ = 0;
  lep1PtErr_     = 999.;
  lep2_       	 = LorentzVector();
  lid2_          = 0;
  lep2IsTightMuon_ = 0;
  lep2IsIsolated_ = 0;
  lep2PtErr_     = 999.;
  lep3_       	 = LorentzVector();
  lid3_          = 0;
  lep3IsTightMuon_ = 0;

  ntaus_          = 0;
  tau1_       	 = LorentzVector();
  tau2_       	 = LorentzVector();

  nphotons_      = 0;
  pho1_       	 = LorentzVector();
  pho2_       	 = LorentzVector();
/*   pho3_       	 = LorentzVector(); */
/*   pho4_       	 = LorentzVector(); */

  njets_ = 0;
  noiseCleaning_  = 0;
  jet1_     = LorentzVector();
  jet1CHF_  = 0.;
  jet1NHF_  = 0.;
  jet1NEMF_  = 0.;
  jet1Unc_  = 0.;
  jet1Btag_ = -999.;
  jet1QGtag_= -999.;  
  jet1PartonId_= 0;  
  jet1QGRho_   = -1.;
  jet1QGPtD_   = -1.;
  jet1QGAxis1_ = -1.;
  jet1QGAxis2_ = -1.;
  jet1QGMult_  = -1.;
  jet2_     = LorentzVector();
  jet2CHF_  = 0.;
  jet2NHF_  = 0.;
  jet2NEMF_  = 0.;
  jet2Unc_  = 0.;
  jet2Btag_ = -999.;
  jet2QGtag_= -999.;
  jet3_     = LorentzVector();
  jet3CHF_  = 0.;
  jet3NHF_  = 0.;
  jet3NEMF_  = 0.;
  jet3Btag_ = -999.;
  jet4_     = LorentzVector();
  jet4Btag_ = -999.;

  ntracks_ = 0;
  track1_ = LorentzVector();
  track2_ = LorentzVector();
  track3_ = LorentzVector();

  Q_		 = -999.;
  id1_  	 = -999.;
  x1_		 = -999.;
  pdf1_ 	 = -999.;  
  id2_  	 = -999.;  
  x2_		 = -999.;
  pdf2_ 	 = -999.;  
  processId_	 = 0;
  npu_           = -999.;
  npuPlusOne_    = -999.;
  npuMinusOne_   = -999.;
  metFiltersWord_= -1.;  
}

#endif
