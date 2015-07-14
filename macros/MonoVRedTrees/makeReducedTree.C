//--------------------------------------------------------------------------------------------------
// Reduce the ntuples and compute all the necessary weights for further studies and plots
//
// Authors: L. Di Matteo                                                                  (Dec 2014)
//--------------------------------------------------------------------------------------------------
#include <iostream>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitMonoJet/Core/interface/MitDMSTree.h"
#include "MitMonoJet/Core/interface/MitLimitTree.h"

#include "TLorentzVector.h"

using namespace std;
using namespace mithep;

//---Cuts Block---//
const double MET_CUT = 180.;
const double MET_CUT_MONOV = 200.;
//-
const double JET_PT_CUT  = 110.;
const double JET_ETA_CUT = 2.5;
const double JET_CHF_CUT = 0.2;  
const double JET_NHF_CUT = 0.7;  
const double JET_NEMF_CUT = 0.7;  
const double JET_NMAX_CUT = 2;
//-
const double BJET_NMAX_CUT = 0;
const double RES_MVA_CUT = 0.6;
//-
const double JET_PT_CUT_INCL = 150.;
const double JET_ETA_CUT_INCL = 2.0;
//-
const double FJET_PT_CUT = 200.;
const double FJET_ETA_CUT = 2.5;
const double T2T1_CUT = 0.5;
const double FJET_MASS_MIN_CUT = 60.;
const double FJET_MASS_MAX_CUT = 110.;
//-
const double JETMET_DPHI_CUT = 2.0;
const double JETJET_DPHI_CUT = 2.0;
//-
const double MU_PT_CUT = 20.;
const double MU_ETA_CUT = 2.1;
const double MU_ID_CUT = 13;
const double MU_TIGHT_CUT = 130;
const double MU_ISO_CUT = 1300;
//-
const double MT_MIN_CUT = 50.;
const double MT_MAX_CUT = 200.;
const double MLL_MIN_CUT = 60.;
const double MLL_MAX_CUT = 120.;
//-
const double PHO_PT_CUT  = 160.;
const double PHO_ETA_CUT = 2.5;

//---Syst block---//
float JESsyst = 0;
bool QGsyst = false;

//---
TString getEnv(const char* name);
//---
void fillOutNtuples(MitLimitTree &outtree, MitDMSTree &intree, double baseWeight, int selMode = 0, bool isData = false, bool exclusive = true, bool testing = false);
//---
bool eventPassSelection(MitDMSTree &intree, int selMode = 0, bool exclusive = true);
//==================================================================================================
void makeReducedTree(int selMode = 0, double lumi = 19700.0, bool updateFile = false, bool exclusive = true, bool testing = false, float setJESsyst = 0, bool setQGsyst = false) 
{
  // Define tree name (depends on selection)
  TString outTreeNameExt;
  if (selMode == 0 || selMode == 4 || selMode == 8) 
    outTreeNameExt = "_signal";
  else if (selMode == 1 || selMode == 5 || selMode == 9) 
    outTreeNameExt = "_di_muon_control";
  else if (selMode == 2 || selMode == 6 || selMode == 10) 
    outTreeNameExt = "_single_muon_control";
  else if (selMode == 3 || selMode == 7 || selMode == 11) 
    outTreeNameExt = "_photon_control";
  else 
    cout << "ERROR -- Incorrect selMode parameter, please review!" << endl;

  // Setup systematics flags
  JESsyst = setJESsyst;
  QGsyst = setQGsyst;

  // Define output files mode
  TString outFileMode = "RECREATE";
  if (updateFile)
    outFileMode = "UPDATE";

  // Read all environment variables
  TString anaDir = getEnv("MIT_MONOJET_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = "boostedv-ana";
  TString prdCfg = getEnv("MIT_PROD_CFG");
  
  // Fix data list for photons
  if (selMode == 3 || selMode == 7 || selMode == 11)
    anaCfg = "boostedv-ana-pj";

  if (testing)
    anaCfg = "boostedv-ana-test";

  // Define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((anaDir + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < samples->NDataSamples(); iSample++) listOfSamples.push_back(samples->GetDataSample(iSample));
  for (UInt_t iSample=0; iSample < samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  

  // Prepare pointer to outfile
  TFile *fin;
  TFile *fout;
  TString outFileName = "boosted";
  if (selMode >= 4 && selMode < 8)
    outFileName = "resolved";
  if (selMode >= 8)
    outFileName = "inclusive";
  if (!exclusive)
    outFileName = "baseline";
  if (testing)
    outFileName = "testing";
    
  if (JESsyst > 0.5)
    outFileName += "_JESup";
  if (JESsyst < -0.5)
    outFileName += "_JESdown";
  
  outFileName += ".root";
  fout = new TFile(outFileName,outFileMode);
  
  // Prepare object to store outtree
  MitLimitTree outTree;
  
  // Generate reduced trees
  // loop through the samples and produce the reduced trees
  for (UInt_t iSample=0; iSample < listOfSamples.size(); iSample++) {
    TString inFilePath = hstDir+"/"+*listOfSamples.at(iSample)->File();
    fin = new TFile(inFilePath,"READ");
    cout << "INFO -- Reading sample file: " << *listOfSamples.at(iSample)->File() << endl;
    // Prepare event weight
    double thisXsec = *listOfSamples.at(iSample)->Xsec();
    double nGenEvts = ((TH1D*)fin->FindObjectAny("hDAllEvents"))->GetEntries();
    double baseWeight = lumi*thisXsec/nGenEvts;
    bool isData = false;
    if (baseWeight < 0)
      isData = true;
    // Read input tree
    MitDMSTree inTree;
    inTree.LoadTree(inFilePath,1);
    inTree.InitTree(1);
    // Start a new sample group according to cfg file legend
    if (*listOfSamples.at(iSample)->Legend() != " ") {
      // Close previous group
      if (iSample > 0) {
        fout->cd();
        outTree.tree_->Write();
        fout->Close();
        fout = TFile::Open(outFileName,"UPDATE");     
      } 
      fout->cd();
      TString outTreeName = *listOfSamples.at(iSample)->Legend();
      outTreeName += outTreeNameExt;
      outTree.CreateTree(outTreeName.Data());
      outTree.InitVariables();
    }
    // Scan on input and fill output ntuple
    cout << "INFO ---> Number of events passing the selection is: ";
    fillOutNtuples(outTree,inTree,baseWeight,selMode,isData,exclusive,testing);
    
    // Close last group
    if (iSample == (listOfSamples.size()-1)) {
      fout->cd();
      outTree.tree_->Write(); 
      fout->Close();
      break;
    }      

    fin->Close();
  }  
  
  return;
}


//==================================================================================================
TString getEnv(const char* name)
{
  if (! gSystem->Getenv(name)) {
    printf(" Environment variable: %s  not defined. EXIT!\n",name);
    return TString("");
  } 
  return TString(gSystem->Getenv(name));  
}

//==================================================================================================
void fillOutNtuples(MitLimitTree &outtree, MitDMSTree &intree, double baseWeight, int selMode, bool isData, bool exclusive, bool testing)
{
  // Prepare special weights if needed
  TFile *fweights;
  TH1F* hweights; 
  TF1* funcweights;
  if (QGsyst) {
    fweights = new TFile("weights.root","READ");
    hweights = (TH1F*) fweights->FindObjectAny("histoWeight");
    funcweights = (TF1*) hweights->GetFunction("fitfunc");
  }

  double weight = -1;
  double sumweight = 0;
  // Loop over tree entries
  Int_t nEntries = intree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
    
    // Get this tree entry
    intree.tree_-> GetEntry(iEntry);

    // Determine if event passes selection
    if (!eventPassSelection(intree,selMode,exclusive))
      continue;
      
    // Determine correctly the event weights
    if (!isData)
      weight = baseWeight*intree.puweight_*intree.genweight_;    
    sumweight += weight;
    
    // Determine the outtree variables for this event, use footprint met for gjets
    outtree.event_ = intree.event_;
    outtree.run_ = intree.run_;
    outtree.lumi_ = intree.lumi_;

    outtree.mvamet_ = intree.met_;
    if (selMode == 1 || selMode == 5 || selMode == 9) {
      outtree.mll_ = intree.mll_;
      outtree.ptll_ = (intree.lep1_ + intree.lep2_).Pt();
    }
    if (selMode == 2 || selMode == 6 || selMode == 10)
      outtree.mt_ = intree.mt_;
    if (selMode == 3 || selMode == 7 || selMode == 11) {
      outtree.mvamet_ = intree.metFprint_;
      outtree.ptpho_ = intree.pho1_.Pt();      
    }
    outtree.mvametphi_ = intree.metPhi_;
    outtree.njets_ = intree.njets_;
    if (JESsyst > 0.5)
      outtree.njets_ = intree.njetsUp_; 
    if (JESsyst < -0.5)
      outtree.njets_ = intree.njetsDown_; 
    outtree.jet1pt_ = intree.fjet1_.Pt()*(1.+JESsyst*intree.fjet1Unc_);
    outtree.genjetpt_ = intree.fjet1_.Pt()*(1.+JESsyst*intree.fjet1Unc_);
    
    if (selMode >= 4 && selMode < 8) {
      outtree.jet1pt_ = (intree.rjet1_ + intree.rjet2_).Pt();
      outtree.genjetpt_ = intree.rjet1_.Pt();
    }
      
    if (selMode >= 8) {
      outtree.jet1pt_ = intree.jet1_.Pt()*(1.+JESsyst*intree.jet1Unc_);
      outtree.genjetpt_ = intree.jet1_.Pt()*(1.+JESsyst*intree.jet1Unc_);
    }
      
    outtree.genVpt_ = intree.genV_.Pt();
    outtree.genVphi_ = intree.genV_.Phi();
    outtree.dmpt_ = intree.genmet_;
    outtree.weight_ = weight;
    if ((selMode == 3 || selMode == 7 || selMode == 11) && intree.genVid_ == 22 && QGsyst)
      outtree.weight_ *= funcweights->Eval(intree.fjet1QGtag_);
    if (isData) {
      outtree.genjetpt_ = -1.;
      outtree.genVpt_ = -1.;
      outtree.genVphi_ = -1.;
      outtree.dmpt_ = -1.;
      outtree.weight_ = -1.;
      if ((selMode == 3 || selMode == 7 || selMode == 11) && QGsyst) 
        outtree.weight_ = funcweights->Eval(intree.fjet1QGtag_);
    }

    // Fill output tree
    outtree.tree_->Fill();
  }
  
  if (isData)
    sumweight = -sumweight;
  if (!isData || testing)
    cout << sumweight << endl;
  else 
    cout << "BLINDED!" << endl;
  return;
}

//==================================================================================================
bool eventPassSelection(MitDMSTree &intree, int selMode, bool exclusive)
{
  // Trigger
  bool triggerBit = (intree.trigger_ & MitDMSTree::HLTJetMet);
  if (selMode == 3 || selMode == 7 || selMode == 11)
    triggerBit = (intree.trigger_ & MitDMSTree::HLTPhoton);

  // Met filters
  bool metFiltersBit = (intree.metFiltersWord_ == 511 || intree.metFiltersWord_ == 1023);
  
  // Preselection
  bool preselBit = (intree.preselWord_ & MitDMSTree::Met); //signal region preselection
  if (selMode == 1 || selMode == 9)
    preselBit = (intree.preselWord_ & MitDMSTree::Zlep); //Z->ll control preselection
  if (selMode == 2 || selMode == 10)
    preselBit = (intree.preselWord_ & MitDMSTree::Wlep); //W->lnu control preselection
  if (selMode == 3 || selMode == 11)
    preselBit = (intree.preselWord_ & MitDMSTree::Gjet); //Photon+jets control preselection
  if (selMode == 4 || selMode == 5 || selMode == 6 || selMode == 7) 
    preselBit = (intree.preselWord_ & MitDMSTree::Resolved); //Resolved preselection

  // Met, use footprint correction for gjets
  float thisMet = intree.met_;
  float thisMetPhi = intree.metPhi_;
  if (selMode == 3 || selMode == 7 || selMode == 11) {
    thisMet = intree.metFprint_;
    thisMetPhi = intree.metFprintPhi_;
  }
  bool metBit = (thisMet > MET_CUT);
  if (selMode < 8)
    metBit = (thisMet > MET_CUT_MONOV);
                      
  // Narrow jets
  float jetpt =  intree.jet1_.Pt()*(1.+JESsyst*intree.jet1Unc_);
  bool jetBit = ((jetpt > JET_PT_CUT && abs(intree.jet1_.eta()) < JET_ETA_CUT)
               && intree.jet1CHF_ > JET_CHF_CUT && intree.jet1NHF_ < JET_NHF_CUT && intree.jet1NEMF_ < JET_NEMF_CUT);
  int njets = intree.njets_;
  if (JESsyst > 0.5)
    njets = intree.njetsUp_; 
  if (JESsyst < -0.5)
    njets = intree.njetsDown_; 
  bool nJetBit = (njets < JET_NMAX_CUT+0.5);

  // Resolved category
  bool resolvedBit = (intree.nbjets_ < BJET_NMAX_CUT+0.5 && intree.rmvaval_ > RES_MVA_CUT);

  // Inclusive category
  bool inclusiveBit = (jetpt > JET_PT_CUT_INCL && abs(intree.jet1_.Eta()) < JET_ETA_CUT_INCL);
  inclusiveBit = inclusiveBit && (abs(MathUtils::DeltaPhi((double)thisMetPhi,intree.jet1_.Phi())) > JETMET_DPHI_CUT);
  inclusiveBit = inclusiveBit && (intree.jet1jet2Dphi_ < -5. || abs(intree.jet1jet2Dphi_) < JETJET_DPHI_CUT);

  // Fat jet category::cut-based
  float fjetpt =  intree.fjet1_.Pt()*(1.+JESsyst*intree.fjet1Unc_);
  bool fatJetBit = ((fjetpt > FJET_PT_CUT && abs(intree.fjet1_.Eta()) < FJET_ETA_CUT)
                  && intree.fjet1CHF_ > JET_CHF_CUT && intree.fjet1NHF_ < JET_NHF_CUT && intree.fjet1NEMF_ < JET_NEMF_CUT);
  fatJetBit = fatJetBit && (abs(MathUtils::DeltaPhi((double)thisMetPhi,intree.fjet1_.Phi())) > JETMET_DPHI_CUT);
  fatJetBit = fatJetBit && (intree.fjet1jet2Dphi_ < -5. || abs(intree.fjet1jet2Dphi_) < JETJET_DPHI_CUT);
  fatJetBit = fatJetBit && (intree.fjet1Tau2_/intree.fjet1Tau1_  < T2T1_CUT 
                            && FJET_MASS_MIN_CUT < intree.fjet1MassPruned_ && intree.fjet1MassPruned_ < FJET_MASS_MAX_CUT
                            && thisMet > MET_CUT_MONOV);

  // Vetoes
  bool vetoBit = (intree.ntaus_ == 0 && intree.nele_ == 0);
  if (selMode == 0 || selMode == 4 || selMode == 8)
    vetoBit = vetoBit && (intree.nphotons_ == 0 && intree.nlep_ == 0);
  if (selMode == 1 || selMode == 5 || selMode == 9)
    vetoBit = vetoBit && (intree.nphotons_ == 0);
  if (selMode == 2 || selMode == 6 || selMode == 10)
    vetoBit = vetoBit && (intree.nphotons_ == 0);
  if (selMode == 3 || selMode == 7 || selMode == 11)
    vetoBit = vetoBit && (intree.nlep_ == 0);
        
  // Extra
  bool extraBit = true;
  //Zll: tight,loose selection
  if (selMode == 1 || selMode == 5 || selMode == 9) {
    extraBit = extraBit && (intree.nlep_ == 2 && ((intree.lid1_*intree.lid2_) <= (MU_ID_CUT*MU_TIGHT_CUT))
                                              && (abs(intree.lep1_.Eta()) < MU_ETA_CUT && abs(intree.lep2_.Eta()) < MU_ETA_CUT));
    extraBit = extraBit && (MLL_MIN_CUT < intree.mll_ && intree.mll_ < MLL_MAX_CUT);    
  } 
  //Wlv
  if (selMode == 2 || selMode == 6 || selMode == 10) {
    extraBit = extraBit && 
              (intree.nlep_ == 1 && abs(intree.lid1_) >= MU_TIGHT_CUT && intree.lep1_.Pt() > MU_PT_CUT && abs(intree.lep1_.Eta()) < MU_ETA_CUT);
    extraBit = extraBit && (intree.mt_ > MT_MIN_CUT && intree.mt_ < MT_MAX_CUT);
  } 
  //Pj
  if (selMode == 3 || selMode == 7 || selMode == 11) {
    extraBit = extraBit && 
              (intree.pho1_.Pt() > PHO_PT_CUT && abs(intree.pho1_.Eta()) < PHO_ETA_CUT);
  } 
        
  //cout 
  //<< triggerBit << " " << metFiltersBit << " " << preselBit << " " << metBit << " " 
  //<< jetBit << " " << nJetBit << " " << resolvedBit << " " << inclusiveBit << " " <<  fatJetBit << " " << vetoBit << " " << extraBit 
  //<< endl;

  bool theDecision;
  // Boosted selection
  if (selMode < 4) 
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && jetBit && nJetBit && fatJetBit && vetoBit && extraBit;
  // Resolved selection and discard boosted events
  else if (selMode >= 4 && selMode < 8)
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && resolvedBit && vetoBit && extraBit && !(fatJetBit&&jetBit&&nJetBit);
  // Inclusive selection and discard boosted/resolved events
  else if (selMode >= 8 && exclusive)
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && jetBit && nJetBit && inclusiveBit && vetoBit && extraBit && !fatJetBit && !resolvedBit;
  // Inclusive selection keep boosted/resolved events
  else if (selMode >= 4 && !exclusive)
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && jetBit && nJetBit && inclusiveBit && vetoBit && extraBit;
  else 
    cout << "ERROR - Selection logic is wrong! Please fix" << endl;
  
  return theDecision;
}
