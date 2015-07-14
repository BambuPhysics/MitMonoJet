#ifndef MitDMSTree_H
#define MitDMSTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

class MitDMSTree {
 public:
  /// float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  /// bit map
  /// DON'T CHANGE ORDER
  enum Trigger {
    HLTMuon     = 1UL<<0,    // event passes single muon trigger
    HLTJetMet   = 1UL<<1,    // event passes jet+met trigger
    HLTVBF      = 1UL<<2,    // event passes VBF trigger
    HLTPhoton   = 1UL<<3     // event passes single photon trigger
  };

  /// bit map
  /// DON'T CHANGE ORDER
  enum HLTMatch {
    JetMatch    = 1UL<<0,    // hardest jet is matched to HLT jet
    MuonMatch   = 1UL<<1,    // hardest lepton is matched to HLT muon
    PhotonMatch = 1UL<<2     // hardest photon is matched to HLT photon
  };

  /// bit map
  /// DON'T CHANGE ORDER
  enum Presel {
    Top      = 1UL<<0,    // event passes top preselection
    Wlep     = 1UL<<1,    // event passes W>lv preselection
    Zlep     = 1UL<<2,    // event passes Z>ll preselection
    Met      = 1UL<<3,    // event passes MET preselection
    Vbf      = 1UL<<4,    // event passes VBF preselection
    Gjet     = 1UL<<5,    // event passes G+jets preselection
    Resolved = 1UL<<6     // event passes Resolved preselection
  };

  /// variables
  unsigned int   isData_;
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   trigger_;
  unsigned int   HLTmatch_;
  unsigned int   nvtx_;
  float          metRaw_;
  float          metRawPhi_;
  float          met_;
  float          metPhi_;
  float          metFprint_;
  float          metFprintPhi_;
  float          genmet_;
  float          genmetPhi_;
  float          mt_;
  float          mll_;

  float          fjet1metDphi_;
  float          fjet1jet2Dphi_;
  float          jet1metDphi_;
  float          jet1jet2Dphi_;
 
  unsigned int   nele_;
  unsigned int   nlep_;
  LorentzVector  lep1_;
  int            lid1_;
  LorentzVector  lep2_;
  int            lid2_;

  int            ntaus_;
  LorentzVector  tau1_;

  unsigned int   nphotons_;
  LorentzVector  pho1_;

  // here come the tag jets (can contain substructure)
  // substructure jets are only saved for the hardest fat jet
  int            nfjets_;
  LorentzVector  fjet1_;
  float          fjet1Unc_;  
  float          fjet1CHF_;  
  float          fjet1NHF_;  
  float          fjet1NEMF_; 
  float          fjet1Btag_;
  float          fjet1Charge_;
  float          fjet1QGtag_;
  float          fjet1Tau1_;
  float          fjet1Tau2_;
  float          fjet1Tau3_;
  float          fjet1C2b0_;
  float          fjet1C2b0p2_;      
  float          fjet1C2b0p5_;      
  float          fjet1C2b1_;        
  float          fjet1C2b2_;        
  float          fjet1QJetVol_;     
  float          fjet1MassSDbm1_;   
  float          fjet1MassSDb0_;    
  float          fjet1MassSDb1_;    
  float          fjet1MassSDb2_;    
  float          fjet1MassPruned_;  
  float          fjet1MassFiltered_;
  float          fjet1MassTrimmed_; 
  float          fjet1Pull_;
  float          fjet1PullAngle_;   
  float          fjet1QGtagSub1_;
  float          fjet1QGPtDSub1_;
  float          fjet1QGAxis1Sub1_;
  float          fjet1QGAxis2Sub1_;
  float          fjet1QGMultSub1_;
  float          fjet1QGtagSub2_;
  float          fjet1QGPtDSub2_;
  float          fjet1QGAxis1Sub2_;
  float          fjet1QGAxis2Sub2_;
  float          fjet1QGMultSub2_;
  unsigned int   fjet1PartonId_;
 
  unsigned int   fjet1nsj_;
  LorentzVector  fjet1sj1_;
  LorentzVector  fjet1sj2_;

  unsigned int   njets_;
  unsigned int   njetsUp_;
  unsigned int   njetsDown_;
  LorentzVector  jet1_;
  float          jet1Unc_;  
  float          jet1CHF_;  
  float          jet1NHF_;  
  float          jet1NEMF_; 
  LorentzVector  jet2_;
  LorentzVector  jet3_;
  LorentzVector  jet4_;
  LorentzVector  jet5_;
 
  unsigned int   nbjets_;
  LorentzVector  bjet1_;
  float          bjet1Btag_;
  LorentzVector  bjet2_;
  float          bjet2Btag_;

  float          rmvaval_;
  LorentzVector  rjet1_;
  float          rjet1_pullang_;
  float          rjet1_qgl_;
  LorentzVector  rjet2_;
  float          rjet2_pullang_;
  float          rjet2_qgl_;
  float          rmdrop_;
  float          rptOverM_;
 
  LorentzVector  genV_;
  unsigned int   genVid_;
  unsigned int   genVdaughterId_;
  float          topPt_;
  float          topBarPt_;
 
  float          genweight_;
  float          Q_;
  float          id1_;
  float          x1_;
  float          pdf1_;
  float          id2_;
  float          x2_;
  float          pdf2_;
  int            processId_;
  float          puweight_;
  float          npu_;
  float          npuPlusOne_;
  float          npuMinusOne_;
  unsigned int   metFiltersWord_;
  unsigned int   preselWord_;

  float          bdt_all_;

 public:
  /// this is the main element
  TFile *f_;
  TTree *tree_;

  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  MitDMSTree():
    lepPtr1_(&lep1_),lepPtr2_(&lep2_),
    tauPtr1_(&tau1_),phoPtr1_(&pho1_),
    fjet1Ptr_(&fjet1_),
    fjet1sjPtr1_(&fjet1sj1_),fjet1sjPtr2_(&fjet1sj2_),
    jetPtr1_(&jet1_),jetPtr2_(&jet2_),jetPtr3_(&jet3_),jetPtr4_(&jet4_),jetPtr5_(&jet5_),
    bjetPtr1_(&bjet1_),bjetPtr2_(&bjet2_),
    rjetPtr1_(&rjet1_),rjetPtr2_(&rjet2_),
    genVPtr_(&genV_) {}
  /// default destructor
  ~MitDMSTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitDMSTree
  void LoadTree(const char* file, int type = -1){
    // to load three different ntuples in the same job DMTree0/1
    // type == 0 standard variables
    // type == 1 bdt_all added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    f_ = TFile::Open(file);
    assert(f_);
    if      (type == 0) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("DMSTree"));
    else if (type == 1) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("DMSTree"));
    else                tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("tree"));
    assert(tree_);
  }

  /// load a MitDMSTree
  void LoadTree(const char* file, const char* treeName){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(treeName));
    assert(tree_);
  }

  /// create a MitDMSTree
  void CreateTree(int type = -1){
    assert(type==type); // just to suppress warnings
    // to create three different ntuples in the same job DMTree0/1
    // type == 0/1 add all variables
    // type = -1 (default) add a minimum set of variables with tree as name
    if     (type == 0) tree_ = new TTree("DMSTree","DM&S ntuples");
    else               tree_ = new TTree("tree","Smurf ntuples");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("isData"          , &isData_          ,   "isData/i");
    tree_->Branch("event"           , &event_           ,   "event/i");
    tree_->Branch("run"             , &run_             ,   "run/i");
    tree_->Branch("lumi"            , &lumi_            ,   "lumi/i");
    tree_->Branch("trigger"         , &trigger_         ,   "trigger/i");
    tree_->Branch("HLTmatch"        , &HLTmatch_        ,   "HLTmatch/i");
    tree_->Branch("nvtx"            , &nvtx_            ,   "nvtx/i");
    tree_->Branch("metRaw"          , &metRaw_          ,   "metRaw/F");
    tree_->Branch("metRawPhi"       , &metRawPhi_       ,   "metRawPhi/F");
    tree_->Branch("met"             , &met_             ,   "met/F");
    tree_->Branch("metPhi"          , &metPhi_          ,   "metPhi/F");
    tree_->Branch("metFprint"       , &metFprint_       ,   "metFprint/F");
    tree_->Branch("metFprintPhi"    , &metFprintPhi_    ,   "metFprintPhi/F");
    tree_->Branch("genmet"          , &genmet_          ,   "genmet/F");
    tree_->Branch("genmetPhi"       , &genmetPhi_       ,   "genmetPhi/F");
    tree_->Branch("mt"              , &mt_              ,   "mt/F");
    tree_->Branch("mll"             , &mll_             ,   "mll/F");

    tree_->Branch("fjet1metDphi"    , &fjet1metDphi_    ,   "fjet1metDphi/F");
    tree_->Branch("fjet1jet2Dphi"   , &fjet1jet2Dphi_   ,   "fjet1jet2Dphi/F");
    tree_->Branch("jet1metDphi"     , &jet1metDphi_     ,   "jet1metDphi/F");
    tree_->Branch("jet1jet2Dphi"    , &jet1jet2Dphi_    ,   "jet1jet2Dphi/F");

    tree_->Branch("nele"         , &nele_         ,   "nele/i");
    tree_->Branch("nlep"         , &nlep_         ,   "nlep/i");
    tree_->Branch("lep1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lid1"         , &lid1_         ,   "lid1/I");
    tree_->Branch("lep2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lid2"         , &lid2_         ,   "lid2/I");

    tree_->Branch("ntaus"        , &ntaus_        ,   "ntaus/i");
    tree_->Branch("tau1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tauPtr1_);

    tree_->Branch("nphotons"     , &nphotons_     ,   "nphotons/i");
    tree_->Branch("pho1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr1_);

    tree_->Branch("nfjets", &nfjets_, "nfjets/i");
    tree_->Branch("fjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &fjet1Ptr_);
    tree_->Branch("fjet1Unc"         , &fjet1Unc_         , "fjet1Unc/F");  
    tree_->Branch("fjet1CHF"         , &fjet1CHF_         , "fjet1CHF/F");  
    tree_->Branch("fjet1NHF"         , &fjet1NHF_         , "fjet1NHF/F");  
    tree_->Branch("fjet1NEMF"        , &fjet1NEMF_        , "fjet1NEMF/F"); 
    tree_->Branch("fjet1Btag"        , &fjet1Btag_        , "fjet1Btag/F");
    tree_->Branch("fjet1Charge"      , &fjet1Charge_      , "fjet1Charge/F");
    tree_->Branch("fjet1QGtag"       , &fjet1QGtag_       , "fjet1QGtag/F");
    tree_->Branch("fjet1Tau1"        , &fjet1Tau1_        , "fjet1Tau1/F");
    tree_->Branch("fjet1Tau2"        , &fjet1Tau2_        , "fjet1Tau2/F");
    tree_->Branch("fjet1Tau3"        , &fjet1Tau3_        , "fjet1Tau3/F");
    tree_->Branch("fjet1C2b0"        , &fjet1C2b0_        , "fjet1C2b0/F");
    tree_->Branch("fjet1C2b0p2"      , &fjet1C2b0p2_      , "fjet1C2b0p2/F");     
    tree_->Branch("fjet1C2b0p5"      , &fjet1C2b0p5_      , "fjet1C2b0p5/F");     
    tree_->Branch("fjet1C2b1"        , &fjet1C2b1_        , "fjet1C2b1/F");       
    tree_->Branch("fjet1C2b2"        , &fjet1C2b2_        , "fjet1C2b2/F");       
    tree_->Branch("fjet1QJetVol"     , &fjet1QJetVol_     , "fjet1QJetVol/F");    
    tree_->Branch("fjet1MassSDbm1"   , &fjet1MassSDbm1_   , "fjet1MassSDbm1/F");  
    tree_->Branch("fjet1MassSDb0"    , &fjet1MassSDb0_    , "fjet1MassSDb0/F");   
    tree_->Branch("fjet1MassSDb1"    , &fjet1MassSDb1_    , "fjet1MassSDb1/F");   
    tree_->Branch("fjet1MassSDb2"    , &fjet1MassSDb2_    , "fjet1MassSDb2/F");   
    tree_->Branch("fjet1MassPruned"  , &fjet1MassPruned_  , "fjet1MassPruned/F"); 
    tree_->Branch("fjet1MassFiltered", &fjet1MassFiltered_, "fjet1MassFiltered/F");
    tree_->Branch("fjet1MassTrimmed" , &fjet1MassTrimmed_ , "fjet1MassTrimmed/F");
    tree_->Branch("fjet1Pull"        , &fjet1Pull_        , "fjet1Pull/F");  
    tree_->Branch("fjet1PullAngle"   , &fjet1PullAngle_   , "fjet1PullAngle/F"); 
    tree_->Branch("fjet1QGtagSub1"   , &fjet1QGtagSub1_   , "fjet1QGtagSub1/F");
    tree_->Branch("fjet1QGPtDSub1"   , &fjet1QGPtDSub1_   , "fjet1QGPtDSub1/F");
    tree_->Branch("fjet1QGAxis1Sub1" , &fjet1QGAxis1Sub1_ , "fjet1QGAxis1Sub1/F");
    tree_->Branch("fjet1QGAxis2Sub1" , &fjet1QGAxis2Sub1_ , "fjet1QGAxis2Sub1/F");
    tree_->Branch("fjet1QGMultSub1"  , &fjet1QGMultSub1_  , "fjet1QGMultSub1/F");
    tree_->Branch("fjet1QGtagSub2"   , &fjet1QGtagSub2_   , "fjet1QGtagSub2/F");
    tree_->Branch("fjet1QGPtDSub2"   , &fjet1QGPtDSub2_   , "fjet1QGPtDSub2/F");
    tree_->Branch("fjet1QGAxis1Sub2" , &fjet1QGAxis1Sub2_ , "fjet1QGAxis1Sub2/F");
    tree_->Branch("fjet1QGAxis2Sub2" , &fjet1QGAxis2Sub2_ , "fjet1QGAxis2Sub2/F");
    tree_->Branch("fjet1QGMultSub2"  , &fjet1QGMultSub2_  , "fjet1QGMultSub2/F");
    tree_->Branch("fjet1PartonId"    , &fjet1PartonId_    , "fjet1PartonId/i");

    tree_->Branch("fjet1nsj", &fjet1nsj_, "fjet1nsj/i");
    tree_->Branch("fjet1sj1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &fjet1sjPtr1_);
    tree_->Branch("fjet1sj2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &fjet1sjPtr2_);

    tree_->Branch("njets"    , &njets_    , "njets/i");
    tree_->Branch("njetsUp"  , &njetsUp_  , "njetsUp/i");
    tree_->Branch("njetsDown", &njetsDown_, "njetsDown/i");
    tree_->Branch("jet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet1Unc"         , &jet1Unc_         , "jet1Unc/F");  
    tree_->Branch("jet1CHF"         , &jet1CHF_         , "jet1CHF/F");  
    tree_->Branch("jet1NHF"         , &jet1NHF_         , "jet1NHF/F");  
    tree_->Branch("jet1NEMF"        , &jet1NEMF_        , "jet1NEMF/F"); 
    tree_->Branch("jet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
    tree_->Branch("jet4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr4_);
    tree_->Branch("jet5", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr5_);

    tree_->Branch("nbjets", &nbjets_, "nbjets/i");
    tree_->Branch("bjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &bjetPtr1_);
    tree_->Branch("bjet1Btag", &bjet1Btag_, "bjet1Btag/F");
    tree_->Branch("bjet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &bjetPtr2_);
    tree_->Branch("bjet2Btag", &bjet2Btag_, "bjet2Btag/F");

    tree_->Branch("rmvaval"      , &rmvaval_      , "rmvaval/F");
    tree_->Branch("rjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &rjetPtr1_);
    tree_->Branch("rjet1_pullang", &rjet1_pullang_, "rjet1_pullang/F");
    tree_->Branch("rjet1_qgl"    , &rjet1_qgl_    , "rjet1_qgl/F");
    tree_->Branch("rjet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &rjetPtr2_);
    tree_->Branch("rjet2_pullang", &rjet2_pullang_, "rjet2_pullang/F");
    tree_->Branch("rjet2_qgl"    , &rjet2_qgl_    , "rjet2_qgl/F");
    tree_->Branch("rmdrop"       , &rmdrop_       , "rmdrop/F");
    tree_->Branch("rptOverM"     , &rptOverM_     , "rptOverM_/F");

    tree_->Branch("genV", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genVPtr_);
    tree_->Branch("genVid",            &genVid_,             "genVid/i");
    tree_->Branch("genVdaughterId",    &genVdaughterId_,     "genVdaughterId/i");
    tree_->Branch("topPt",             &topPt_    ,          "topPt/F");
    tree_->Branch("topBarPt",          &topBarPt_ ,          "topBarPt/F");

    tree_->Branch("genweight",      &genweight_   ,   "genweight/F");
    tree_->Branch("Q",              &Q_,              "Q/F");
    tree_->Branch("id1",            &id1_,            "id1/F");
    tree_->Branch("x1",             &x1_,             "x1/F");
    tree_->Branch("pdf1",           &pdf1_,           "pdf1/F");
    tree_->Branch("id2",            &id2_,            "id2/F");
    tree_->Branch("x2",             &x2_,             "x2/F");
    tree_->Branch("pdf2",           &pdf2_,           "pdf2/F");
    tree_->Branch("processId",      &processId_ ,     "processId/I");
    tree_->Branch("puweight",       &puweight_,       "puweight/F");
    tree_->Branch("npu",            &npu_,            "npu/F");
    tree_->Branch("npuPlusOne",     &npuPlusOne_,     "npuPlusOne/F");
    tree_->Branch("npuMinusOne",    &npuMinusOne_,    "npuMinusOne/F");
    tree_->Branch("metFiltersWord", &metFiltersWord_, "metFiltersWord/I");
    tree_->Branch("preselWord",     &preselWord_,     "preselWord/I");

  }

  // initialze a MitDMSTree
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
    tree_->SetBranchAddress("isData",        &isData_);
    tree_->SetBranchAddress("event",         &event_);
    tree_->SetBranchAddress("run",           &run_);
    tree_->SetBranchAddress("lumi",          &lumi_);
    tree_->SetBranchAddress("trigger",       &trigger_);
    tree_->SetBranchAddress("HLTmatch",      &HLTmatch_);
    tree_->SetBranchAddress("nvtx",          &nvtx_);
    tree_->SetBranchAddress("metRaw",        &metRaw_);
    tree_->SetBranchAddress("metRawPhi",     &metRawPhi_);
    tree_->SetBranchAddress("met",           &met_);
    tree_->SetBranchAddress("metPhi",        &metPhi_);
    tree_->SetBranchAddress("metFprint",     &metFprint_);
    tree_->SetBranchAddress("metFprintPhi",  &metFprintPhi_);
    tree_->SetBranchAddress("genmet",        &genmet_);
    tree_->SetBranchAddress("genmetPhi",     &genmetPhi_);
    tree_->SetBranchAddress("mt",            &mt_);
    tree_->SetBranchAddress("mll",           &mll_);

    tree_->SetBranchAddress("fjet1metDphi" , &fjet1metDphi_);
    tree_->SetBranchAddress("fjet1jet2Dphi", &fjet1jet2Dphi_);
    tree_->SetBranchAddress("jet1metDphi"  , &jet1metDphi_);
    tree_->SetBranchAddress("jet1jet2Dphi" , &jet1jet2Dphi_);

    tree_->SetBranchAddress("nele",          &nele_);
    tree_->SetBranchAddress("nlep",          &nlep_);
    tree_->SetBranchAddress("lep1",          &lepPtr1_);
    tree_->SetBranchAddress("lid1",          &lid1_);
    tree_->SetBranchAddress("lep2",          &lepPtr2_);
    tree_->SetBranchAddress("lid2",          &lid2_);

    tree_->SetBranchAddress("ntaus",         &ntaus_);
    tree_->SetBranchAddress("tau1",          &tauPtr1_);

    tree_->SetBranchAddress("nphotons",      &nphotons_);
    tree_->SetBranchAddress("pho1",          &phoPtr1_);

    tree_->SetBranchAddress("nfjets"           , &nfjets_           );
    tree_->SetBranchAddress("fjet1"            , &fjet1Ptr_         );
    tree_->SetBranchAddress("fjet1Unc"         , &fjet1Unc_         );  
    tree_->SetBranchAddress("fjet1CHF"         , &fjet1CHF_          );  
    tree_->SetBranchAddress("fjet1NHF"         , &fjet1NHF_          );  
    tree_->SetBranchAddress("fjet1NEMF"        , &fjet1NEMF_         ); 
    tree_->SetBranchAddress("fjet1Btag"        , &fjet1Btag_        );
    tree_->SetBranchAddress("fjet1Charge"      , &fjet1Charge_      );
    tree_->SetBranchAddress("fjet1QGtag"       , &fjet1QGtag_       );
    tree_->SetBranchAddress("fjet1Tau1"        , &fjet1Tau1_        );
    tree_->SetBranchAddress("fjet1Tau2"        , &fjet1Tau2_        );
    tree_->SetBranchAddress("fjet1Tau3"        , &fjet1Tau3_        );
    tree_->SetBranchAddress("fjet1C2b0"        , &fjet1C2b0_        );
    tree_->SetBranchAddress("fjet1C2b0p2"      , &fjet1C2b0p2_      );     
    tree_->SetBranchAddress("fjet1C2b0p5"      , &fjet1C2b0p5_      );     
    tree_->SetBranchAddress("fjet1C2b1"        , &fjet1C2b1_        );       
    tree_->SetBranchAddress("fjet1C2b2"        , &fjet1C2b2_        );       
    tree_->SetBranchAddress("fjet1QJetVol"     , &fjet1QJetVol_     );    
    tree_->SetBranchAddress("fjet1MassSDbm1"   , &fjet1MassSDbm1_   );  
    tree_->SetBranchAddress("fjet1MassSDb0"    , &fjet1MassSDb0_    );   
    tree_->SetBranchAddress("fjet1MassSDb1"    , &fjet1MassSDb1_    );   
    tree_->SetBranchAddress("fjet1MassSDb2"    , &fjet1MassSDb2_    );   
    tree_->SetBranchAddress("fjet1MassPruned"  , &fjet1MassPruned_  ); 
    tree_->SetBranchAddress("fjet1MassFiltered", &fjet1MassFiltered_);
    tree_->SetBranchAddress("fjet1MassTrimmed" , &fjet1MassTrimmed_ );
    tree_->SetBranchAddress("fjet1Pull"        , &fjet1Pull_        );  
    tree_->SetBranchAddress("fjet1PullAngle"   , &fjet1PullAngle_   ); 
    tree_->SetBranchAddress("fjet1QGtagSub1"   , &fjet1QGtagSub1_   );
    tree_->SetBranchAddress("fjet1QGPtDSub1"   , &fjet1QGPtDSub1_   );
    tree_->SetBranchAddress("fjet1QGAxis1Sub1" , &fjet1QGAxis1Sub1_ );
    tree_->SetBranchAddress("fjet1QGAxis2Sub1" , &fjet1QGAxis2Sub1_ );
    tree_->SetBranchAddress("fjet1QGMultSub1"  , &fjet1QGMultSub1_  );
    tree_->SetBranchAddress("fjet1QGtagSub2"   , &fjet1QGtagSub2_   );
    tree_->SetBranchAddress("fjet1QGPtDSub2"   , &fjet1QGPtDSub2_   );
    tree_->SetBranchAddress("fjet1QGAxis1Sub2" , &fjet1QGAxis1Sub2_ );
    tree_->SetBranchAddress("fjet1QGAxis2Sub2" , &fjet1QGAxis2Sub2_ );
    tree_->SetBranchAddress("fjet1QGMultSub2"  , &fjet1QGMultSub2_  );
    tree_->SetBranchAddress("fjet1PartonId"    , &fjet1PartonId_    );

    tree_->SetBranchAddress("fjet1nsj"         , &fjet1nsj_         );
    tree_->SetBranchAddress("fjet1sj1"         , &fjet1sjPtr1_      );
    tree_->SetBranchAddress("fjet1sj2"         , &fjet1sjPtr2_      );

    tree_->SetBranchAddress("njets"            , &njets_            );
    tree_->SetBranchAddress("njetsUp"          , &njetsUp_          );
    tree_->SetBranchAddress("njetsDown"        , &njetsDown_        );
    tree_->SetBranchAddress("jet1"             , &jetPtr1_          );
    tree_->SetBranchAddress("jet1Unc"          , &jet1Unc_          );  
    tree_->SetBranchAddress("jet1CHF"          , &jet1CHF_          );  
    tree_->SetBranchAddress("jet1NHF"          , &jet1NHF_          );  
    tree_->SetBranchAddress("jet1NEMF"         , &jet1NEMF_         ); 
    tree_->SetBranchAddress("jet2"             , &jetPtr2_          );
    tree_->SetBranchAddress("jet3"             , &jetPtr3_          );
    tree_->SetBranchAddress("jet4"             , &jetPtr4_          );
    tree_->SetBranchAddress("jet5"             , &jetPtr5_          );

    tree_->SetBranchAddress("nbjets"           , &nbjets_           );
    tree_->SetBranchAddress("bjet1"            , &bjetPtr1_         );
    tree_->SetBranchAddress("bjet1Btag"        , &bjet1Btag_        );
    tree_->SetBranchAddress("bjet2"            , &bjetPtr2_         );
    tree_->SetBranchAddress("bjet2Btag"        , &bjet2Btag_        );

    tree_->SetBranchAddress("rmvaval"          , &rmvaval_          );
    tree_->SetBranchAddress("rjet1"            , &rjetPtr1_         );
    tree_->SetBranchAddress("rjet1_pullang"    , &rjet1_pullang_    );
    tree_->SetBranchAddress("rjet1_qgl"        , &rjet1_qgl_        );
    tree_->SetBranchAddress("rjet2"            , &rjetPtr2_         );
    tree_->SetBranchAddress("rjet2_pullang"    , &rjet2_pullang_    );
    tree_->SetBranchAddress("rjet2_qgl"        , &rjet2_qgl_        );
    tree_->SetBranchAddress("rmdrop"           , &rmdrop_           );
    tree_->SetBranchAddress("rptOverM"         , &rptOverM_         );

    tree_->SetBranchAddress("genV"             , &genVPtr_          );
    tree_->SetBranchAddress("genVid "          , &genVid_           );
    tree_->SetBranchAddress("genVdaughterId"   , &genVdaughterId_   );
    tree_->SetBranchAddress("topPt"            , &topPt_            );
    tree_->SetBranchAddress("topBarPt"         , &topBarPt_         );

    tree_->SetBranchAddress("Q"             ,	&Q_             );
    tree_->SetBranchAddress("id1"           ,	&id1_           );
    tree_->SetBranchAddress("x1"            ,	&x1_            );
    tree_->SetBranchAddress("pdf1"          ,	&pdf1_          );
    tree_->SetBranchAddress("id2"           ,	&id2_           );
    tree_->SetBranchAddress("x2"            ,	&x2_            );
    tree_->SetBranchAddress("pdf2"          ,	&pdf2_          );
    tree_->SetBranchAddress("processId"     , &processId_     );
    tree_->SetBranchAddress("puweight"      , &puweight_      );
    tree_->SetBranchAddress("npu"           ,	&npu_           );
    tree_->SetBranchAddress("npuPlusOne"    , &npuPlusOne_    );
    tree_->SetBranchAddress("npuMinusOne"   , &npuMinusOne_   );
    tree_->SetBranchAddress("metFiltersWord", &metFiltersWord_);
    tree_->SetBranchAddress("preselWord"    , &preselWord_    );

    tree_->SetBranchAddress("genweight",      &genweight_     );
    tree_->SetBranchAddress("Q",              &Q_             );
    tree_->SetBranchAddress("id1",            &id1_           );
    tree_->SetBranchAddress("x1",             &x1_            );
    tree_->SetBranchAddress("pdf1",           &pdf1_          );
    tree_->SetBranchAddress("id2",            &id2_           );
    tree_->SetBranchAddress("x2",             &x2_            );
    tree_->SetBranchAddress("pdf2",           &pdf2_          );
    tree_->SetBranchAddress("processId",      &processId_     );
    tree_->SetBranchAddress("puweight",       &puweight_      );
    tree_->SetBranchAddress("npu",            &npu_           );
    tree_->SetBranchAddress("npuPlusOne",     &npuPlusOne_    );
    tree_->SetBranchAddress("npuMinusOne",    &npuMinusOne_   );
    tree_->SetBranchAddress("metFiltersWord", &metFiltersWord_);
    tree_->SetBranchAddress("preselWord",     &preselWord_    );

    if (type == 1) tree_->SetBranchAddress("bdt_all"    , &bdt_all_);

    gErrorIgnoreLevel = currentState;
  }

  private:

  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* tauPtr1_;
  LorentzVector* phoPtr1_;

  LorentzVector* fjet1Ptr_;
  LorentzVector* fjet1sjPtr1_;
  LorentzVector* fjet1sjPtr2_;

  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* jetPtr3_;
  LorentzVector* jetPtr4_;
  LorentzVector* jetPtr5_;

  LorentzVector* bjetPtr1_;
  LorentzVector* bjetPtr2_;

  LorentzVector* rjetPtr1_;
  LorentzVector* rjetPtr2_;

  LorentzVector* genVPtr_;
}; 

inline void 
MitDMSTree::InitVariables(){
  // inizialize variables
  isData_        = 0;
  event_         = 0;
  run_           = 0;
  lumi_          = 0;
  trigger_       = 0;
  HLTmatch_      = 0;
  nvtx_          = 0;
  metRaw_        = -1.;
  metRawPhi_     = -10.;
  met_           = -1.;
  metPhi_        = -10.;
  metFprint_     = -1.;
  metFprintPhi_  = -10.;
  genmet_        = -1.;
  genmetPhi_     = -10.;
  mt_            = -1.;
  mll_           = -1.;

  fjet1metDphi_  = -10.;
  fjet1jet2Dphi_ = -10.;
  jet1metDphi_   = -10.;
  jet1jet2Dphi_  = -10.;

  nele_          = 0;
  nlep_          = 0;
  lep1_       	 = LorentzVector();
  lid1_          = 0;
  lep2_       	 = LorentzVector();
  lid2_          = 0;

  ntaus_         = 0;
  tau1_       	 = LorentzVector();
  
  nphotons_      = 0;
  pho1_       	 = LorentzVector();

  nfjets_         = 0;
  fjet1_          = LorentzVector();
  fjet1Unc_       = -1.;
  fjet1CHF_       = -1.;
  fjet1NHF_       = -1.;
  fjet1NEMF_      = -1.;
  fjet1Btag_      = -1.;
  fjet1Charge_    = -10.;
  fjet1QGtag_     = -1.;
  fjet1Tau1_      = -1.;
  fjet1Tau2_      = -1.;
  fjet1Tau3_      = -1.;
  fjet1C2b0_      = -1.;
  fjet1C2b0p2_      = -1.;
  fjet1C2b0p5_      = -1.;
  fjet1C2b1_        = -1.;
  fjet1C2b2_        = -1.;
  fjet1QJetVol_     = -1.;
  fjet1MassSDbm1_   = -1.;
  fjet1MassSDb0_    = -1.;
  fjet1MassSDb1_    = -1.;
  fjet1MassSDb2_    = -1.;
  fjet1MassPruned_  = -1.;
  fjet1MassFiltered_= -1.;
  fjet1MassTrimmed_ = -1.;
  fjet1Pull_        = -1.;
  fjet1PullAngle_   = -10.;
  fjet1QGtagSub1_   = -1.;
  fjet1QGPtDSub1_   = -1.;
  fjet1QGAxis1Sub1_ = -1.;
  fjet1QGAxis2Sub1_ = -1.;
  fjet1QGMultSub1_  = -1.;
  fjet1QGtagSub2_   = -1.;
  fjet1QGPtDSub2_   = -1.;
  fjet1QGAxis1Sub2_ = -1.;
  fjet1QGAxis2Sub2_ = -1.;
  fjet1QGMultSub2_  = -1.;
  fjet1PartonId_  = 0;

  fjet1nsj_      = 0;
  fjet1sj1_      = LorentzVector();
  fjet1sj2_      = LorentzVector();
  
  njets_         = 0;
  njetsUp_       = 0;
  njetsDown_     = 0;
  jet1_          = LorentzVector();
  jet1Unc_       = -1.;
  jet1CHF_       = -1.;
  jet1NHF_       = -1.;
  jet1NEMF_      = -1.;
  jet2_          = LorentzVector();
  jet3_          = LorentzVector();
  jet4_          = LorentzVector();
  jet5_          = LorentzVector();

  nbjets_        = 0; 
  bjet1_         = LorentzVector();
  bjet1Btag_     = 0;
  bjet2_         = LorentzVector();
  bjet2Btag_     = 0;
  
  rmvaval_      = -10.; 
  rjet1_        = LorentzVector();
  rjet1_pullang_= -10.;
  rjet1_qgl_    = -1.;
  rjet2_        = LorentzVector();
  rjet2_pullang_= -10.;
  rjet2_qgl_    = -1.;
  rmdrop_       = -1.;
  rptOverM_     = -1.;

  genV_          = LorentzVector();
  genVid_        = 0;
  genVdaughterId_= 0;
  topPt_         = 0;
  topBarPt_      = 0;
  
  genweight_     = 1.;
  Q_		         = -1.;
  id1_  	       = -1.;
  x1_		         = -1.;
  pdf1_ 	       = -1.;  
  id2_  	       = -1.;  
  x2_		         = -1.;
  pdf2_ 	       = -1.;  
  processId_	   = 0;
  puweight_      = 1.;
  npu_           = -1.;
  npuPlusOne_    = -1.;
  npuMinusOne_   = -1.;
  metFiltersWord_= 0;  
  preselWord_    = 0;  
  
  bdt_all_ = -10.;
}

#endif
