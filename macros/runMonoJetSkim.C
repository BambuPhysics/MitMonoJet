#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"
#include "MitPhysics/Mods/interface/MuonIDModRun1.h"
#include "MitPhysics/Mods/interface/ElectronIDModRun1.h"
#include "MitPhysics/Mods/interface/PhotonIDModRun1.h"
#include "MitPhysics/Mods/interface/PFTauIDModRun1.h"
#include "MitPhysics/Mods/interface/JetIdMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitMonoJet/SelMods/interface/MonoJetAnalysisMod.h"
//#include "MitTest/Mods/interface/EventCounterMod.h"
#endif

using namespace mithep;

//--------------------------------------------------------------------------------------------------
void runMonoJetSkim(const char *fileset    = "0000",
		    const char *skim       = "noskim",
		    const char *dataset    = "s12-zjets-ptz100-v7a",
		    const char *book       = "t2mit/filefi/032",
		    const char *catalogDir = "/home/cmsprod/catalog",
		    const char *outputLabel = "monojet",
		    int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // json parameters get passed through the environment
  // for MC, the value must be "~"
  //------------------------------------------------------------------------------------------------
  TString json(gSystem->Getenv("MIT_PROD_JSON"));
  if (json.Length() == 0) {
    printf(" JSON file was not properly defined. EXIT!\n");
    return;
  }

  TString jsonFile = TString("/home/cmsprod/cms/json/") + json;
  Bool_t  isData   = (json != "~");

  TString MitData(gSystem->Getenv("MIT_DATA"));
  if(MitData.Length() == 0){
    MitData = gSystem->Getenv("CMSSW_BASE");
    MitData += "/src/MitPhysics/data";
  }

  printf("\n Initialization worked: \n\n");
  printf("   JSON   : %s (file: %s)\n",  json.Data(), jsonFile.Data());
  printf("   isData : %d\n\n",isData);

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //  mithep::EventCounterMod *counter = new mithep::EventCounterMod;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/cmsprod/cms/json/-") != 0)   ) {
    printf("\n Jason file added: %s \n\n",jsonFile.Data());
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }
  printf("\n Run lumi worked. \n\n");

  // Generator info
  GeneratorMod *generatorMod = new GeneratorMod;
  generatorMod->SetPrintDebug(kFALSE);
  generatorMod->SetPtLeptonMin(0.0);
  generatorMod->SetEtaLeptonMax(2.7);
  generatorMod->SetPtPhotonMin(0.0);
  generatorMod->SetEtaPhotonMax(2.7);
  generatorMod->SetPtRadPhotonMin(0.0);
  generatorMod->SetEtaRadPhotonMax(2.7);
  generatorMod->SetIsData(isData);
  generatorMod->SetFillHist(! isData);
  generatorMod->SetApplyISRFilter(kFALSE);
  generatorMod->SetApplyVVFilter(kFALSE);
  generatorMod->SetApplyVGFilter(kFALSE);
  generatorMod->SetFilterBTEvents(kFALSE);

  //-----------------------------------------------------------------------------------------------------------
  // HLT information : trigger not applied (neither for data nor for MC, store info to apply selection offline
  //-----------------------------------------------------------------------------------------------------------
  HLTMod *hltModP = new HLTMod("HLTModP");

  // monojet triggers
  const int nMjtTrigs = 12;
  TString monoJetTriggers[nMjtTrigs] = { "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4",
					 "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3",
					 "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1",
					 "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5",
					 "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4",
					 "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3",
					 "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2",
					 "HLT_MET120_HBHENoiseCleaned_v6",
					 "HLT_MET120_HBHENoiseCleaned_v5",
					 "HLT_MET120_HBHENoiseCleaned_v4",
					 "HLT_MET120_HBHENoiseCleaned_v3",
					 "HLT_MET120_HBHENoiseCleaned_v2" };

  for (int i=0; i<nMjtTrigs; i++)
    hltModP->AddTrigger(TString("!+"+monoJetTriggers[i]),0,999999);

  // VBF triggers
  const int nVbfTrigs = 7;
  TString vbfTriggers[nVbfTrigs] = { "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9",
				     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v8",
				     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v6",
				     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v5",
				     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v4",
				     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v3",
				     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v2" };

  for (int i=0; i<nVbfTrigs; i++)
    hltModP->AddTrigger((TString("!+")+vbfTriggers[i]).Data(),0,999999);

  hltModP->SetBitsName("HLTBits");
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
  hltModP->SetPrintTable(kFALSE);

  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kFALSE);
  
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPvMod = new GoodPVFilterMod;
  goodPvMod->SetMinVertexNTracks(0);
  goodPvMod->SetMinNDof(4.0);
  goodPvMod->SetMaxAbsZ(24.0);
  goodPvMod->SetMaxRho(2.0);
  goodPvMod->SetIsMC(!isData);
  goodPvMod->SetVertexesName("PrimaryVertexes");
  
  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------

  //-----------------------------------
  // Lepton Selection 
  //-----------------------------------
  ElectronIDModRun1* eleIdMod = new ElectronIDModRun1;
  eleIdMod->SetPtMin(10.);  
  eleIdMod->SetEtaMax(2.5);
  eleIdMod->SetApplyEcalFiducial(true);
  //eleIdMod->SetIdType(mithep::ElectronTools::kVBTFWorkingPoint95Id);
  eleIdMod->SetIDType("VBTFWorkingPoint95Id");
  eleIdMod->SetIsoType("PFIso");
  eleIdMod->SetApplyConversionFilterType1(kTRUE);
  eleIdMod->SetApplyConversionFilterType2(kFALSE);
  eleIdMod->SetChargeFilter(kFALSE);
  eleIdMod->SetApplyD0Cut(kTRUE);
  eleIdMod->SetApplyDZCut(kTRUE);
  eleIdMod->SetWhichVertex(-1);
  //eleIdMod->SetInputName("Electrons");
  //eleIdMod->SetOutputName("GoodElectronsBS");
  eleIdMod->SetGoodElectronsName("GoodElectronsBS");
  //eleIdMod->SetRhoAlgo(mithep::PileupEnergyDensity::kKt6PFJets);

  MuonIDModRun1 *muonId = new MuonIDModRun1;
  muonId->SetOutputName("GoodMuons");
  muonId->SetIntRadius(0.0);
  muonId->SetClassType("GlobalTracker");
  muonId->SetIDType("WWMuIdV4");
  muonId->SetIsoType("IsoRingsV0_BDTG_Iso");
  muonId->SetApplyD0Cut(kTRUE);
  muonId->SetApplyDZCut(kTRUE);
  muonId->SetWhichVertex(0);
  muonId->SetRhoAlgo(mithep::PileupEnergyDensity::kKt6PFJets);
  muonId->SetPtMin(10.);
  muonId->SetEtaCut(2.4);

  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  electronCleaning->SetCleanMuonsName(muonId->GetOutputName());
  electronCleaning->SetGoodElectronsName(eleIdMod->GetOutputName());
  electronCleaning->SetCleanElectronsName("CleanElectrons");

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName(muonId->GetOutputName());
  merger->SetElectronsName(electronCleaning->GetOutputName());
  merger->SetMergedName("MergedLeptons");

  //-----------------------------------
  // Photon Regression + ID 
  //-----------------------------------
  PhotonMvaMod *photonReg = new PhotonMvaMod;
  photonReg->SetRegressionVersion(3);
  photonReg->SetRegressionWeights((MitData+TString("/gbrv3ph_52x.root")).Data());
  photonReg->SetOutputName("GoodPhotonsRegr");
  photonReg->SetApplyShowerRescaling(kTRUE);
  photonReg->SetMinNumPhotons(0);
  photonReg->SetIsData(isData);

  PhotonIDModRun1 *photonIDMod = new PhotonIDModRun1;
  photonIDMod->SetPtMin(0.0);
  photonIDMod->SetOutputName("GoodPhotons");
  photonIDMod->SetIDType("BaseLineCiCPFNoPresel");
  photonIDMod->SetIsoType("NoIso");
  photonIDMod->SetApplyElectronVeto(kTRUE);
  photonIDMod->SetApplyPixelSeed(kTRUE);
  photonIDMod->SetApplyConversionId(kTRUE);
  photonIDMod->SetApplyFiduciality(kTRUE);       
  photonIDMod->SetIsData(isData);
  photonIDMod->SetPhotonsFromBranch(kFALSE);
  photonIDMod->SetInputName(photonReg->GetOutputName());
  //get the photon with regression energy  
  photonIDMod->DoMCSmear(kTRUE);
  photonIDMod->DoDataEneCorr(kTRUE);
  //---------------------------------shower shape scale--------------------------------------------------------------------------------
  photonIDMod->SetDoShowerShapeScaling(kTRUE);
  photonIDMod->SetShowerShapeType("2012ShowerShape");

  PFTauIDModRun1 *pfTauIDMod = new PFTauIDModRun1;
  pfTauIDMod->SetPFTausName("HPSTaus");
  pfTauIDMod->SetIsLooseId(kFALSE);

  PhotonCleaningMod *photonCleaningMod = new PhotonCleaningMod;
  photonCleaningMod->SetCleanElectronsName(electronCleaning->GetOutputName());
  photonCleaningMod->SetGoodPhotonsName(photonIDMod->GetOutputName());
  photonCleaningMod->SetCleanPhotonsName("CleanPhotons");

  PFTauCleaningMod *pfTauCleaningMod = new PFTauCleaningMod;
  pfTauCleaningMod->SetGoodPFTausName(pfTauIDMod->GetGoodPFTausName());
  pfTauCleaningMod->SetCleanMuonsName(muonId->GetOutputName());

  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if (isData){ 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data());
  }                                                                                      
  else {                                                                                 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    

  JetIdMod *jetID = new JetIdMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetOutputName("GoodJets");
  jetID->SetPtMin(30.0);
  jetID->SetEtaMax(4.7);
  jetID->SetApplyPFLooseId(kTRUE);
  jetID->SetApplyBetaCut(kFALSE);
  jetID->SetMVATrainingSet(JetIDMVA::k53);
  jetID->SetMVACutWP(JetIDMVA::kLoose);
  jetID->SetMVACutsFile(MitData + "/jetIDCuts_121221.dat");
  jetID->SetMVAWeightsFile(MitData + "/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml");

  JetCleaningMod *jetCleaning = new JetCleaningMod;
  jetCleaning->SetCleanElectronsName(electronCleaning->GetOutputName());
  jetCleaning->SetCleanMuonsName(muonId->GetOutputName());
  jetCleaning->SetCleanPhotonsName(photonCleaningMod->GetOutputName());
  jetCleaning->SetApplyPhotonRemoval(kTRUE);
  jetCleaning->SetGoodJetsName(jetID->GetOutputName());
  jetCleaning->SetCleanJetsName("CleanJets");
        
  MetCorrectionMod *metCorrT0T1Shift = new MetCorrectionMod;
  metCorrT0T1Shift->SetInputName("PFMet");
  metCorrT0T1Shift->SetJetsName(pubJet->GetOutputName());    
  metCorrT0T1Shift->SetCorrectedJetsName(jetCorr->GetOutputName());    
  metCorrT0T1Shift->SetCorrectedName("PFMetT0T1Shift");   
  metCorrT0T1Shift->ApplyType0(kTRUE);   
  metCorrT0T1Shift->ApplyType1(kTRUE);   
  metCorrT0T1Shift->ApplyShift(kTRUE);   
  metCorrT0T1Shift->IsData(isData);
  metCorrT0T1Shift->SetPrint(kFALSE);

  //------------------------------------------------------------------------------------------------
  // select events
  //------------------------------------------------------------------------------------------------

  float minLeadingJetEt = 100.;
  float maxJetEta       = 4.7;
  float minMet          = 160.;


  MonoJetAnalysisMod *monojetSel = new MonoJetAnalysisMod("MonoJetSelector");
  monojetSel->SetInputMetName(metCorrT0T1Shift->GetOutputName());
  monojetSel->SetMetFromBranch(kFALSE);
  monojetSel->SetJetsName(jetCleaning->GetOutputName());
  monojetSel->SetJetsFromBranch(kFALSE);
  monojetSel->SetElectronsName(electronCleaning->GetOutputName());
  monojetSel->SetElectronsFromBranch(kFALSE);
  monojetSel->SetMuonsName(muonId->GetOutputName());
  monojetSel->SetMuonsFromBranch(kFALSE);
  monojetSel->SetTausName(pfTauCleaningMod->GetOutputName());
  monojetSel->SetTausFromBranch(kFALSE);
  monojetSel->SetLeptonsName(merger->GetOutputName());
  monojetSel->SetCategoriesName("MonoJetEventCategories");

 
  // Jet + MET (signal region)
  unsigned iCat = 0;
  monojetSel->SetMinNumLeptons(iCat, 0);
  monojetSel->SetMaxNumLeptons(iCat, 0);
  monojetSel->SetMinNumTaus(iCat, 0);
  monojetSel->SetMaxNumTaus(iCat, 0);
  monojetSel->SetMinNumJets(iCat, 1);
  monojetSel->SetMaxNumJets(iCat, 1);
  monojetSel->SetMinNumGenNeutrinos(iCat, 0);
  monojetSel->SetMinJetEt(iCat, minLeadingJetEt);
  monojetSel->SetMaxJetEta(iCat, maxJetEta);
  monojetSel->SetMinMetEt(iCat, minMet);
  monojetSel->SetMinEmulMetEt(iCat, 0.);
  monojetSel->SetMinChargedHadronFrac(iCat, 0.2); 
  monojetSel->SetMaxNeutralHadronFrac(iCat, 0.7);
  monojetSel->SetMaxNeutralEmFrac(iCat, 0.7);

  // Dilepton (Z->ll)
  iCat = 1;
  monojetSel->SetMinNumLeptons(iCat, 2);
  monojetSel->SetMaxNumLeptons(iCat, 2);
  monojetSel->SetMinNumTaus(iCat, 0);
  monojetSel->SetMaxNumTaus(iCat, 0);
  monojetSel->SetMinNumJets(iCat, 1);
  monojetSel->SetMaxNumJets(iCat, 1);
  monojetSel->SetMinNumGenNeutrinos(iCat, 0);
  monojetSel->SetMinJetEt(iCat, minLeadingJetEt);
  monojetSel->SetMaxJetEta(iCat, maxJetEta);
  monojetSel->SetMinMetEt(iCat, 0.);
  monojetSel->SetMinEmulMetEt(iCat, minMet);
  monojetSel->SetMinChargedHadronFrac(iCat, 0.2); 
  monojetSel->SetMaxNeutralHadronFrac(iCat, 0.7);
  monojetSel->SetMaxNeutralEmFrac(iCat, 0.7);
  
  // Single lepton (W->lnu)
  iCat = 2;
  monojetSel->SetMinNumLeptons(iCat, 1);
  monojetSel->SetMaxNumLeptons(iCat, 1);
  monojetSel->SetMinNumTaus(iCat, 0);
  monojetSel->SetMaxNumTaus(iCat, 0);
  monojetSel->SetMinNumJets(iCat, 1);
  monojetSel->SetMaxNumJets(iCat, 1);
  monojetSel->SetMinNumGenNeutrinos(iCat, 0);
  monojetSel->SetMinJetEt(iCat, minLeadingJetEt);
  monojetSel->SetMaxJetEta(iCat, maxJetEta);
  monojetSel->SetMinMetEt(iCat, minMet);
  monojetSel->SetMinEmulMetEt(iCat, 0.);
  monojetSel->SetMinChargedHadronFrac(iCat, 0.2); 
  monojetSel->SetMaxNeutralHadronFrac(iCat, 0.7);
  monojetSel->SetMaxNeutralEmFrac(iCat, 0.7);

  // if dataset is Z MC, skim based on the existence of neutrinos
  if(TString(dataset).Contains("zjets")){
    iCat = 3;
    monojetSel->SetMinNumLeptons(iCat, 0);
    monojetSel->SetMaxNumLeptons(iCat, 0);
    monojetSel->SetMinNumTaus(iCat, 0);
    monojetSel->SetMaxNumTaus(iCat, 0);
    monojetSel->SetMinNumJets(iCat, 1);
    monojetSel->SetMaxNumJets(iCat, 1);
    // this is what makes this module special
    monojetSel->SetMinNumGenNeutrinos(iCat, 2);
    monojetSel->SetMinJetEt(iCat, minLeadingJetEt);
    monojetSel->SetMaxJetEta(iCat, maxJetEta);
    monojetSel->SetMinMetEt(iCat, minMet);
    monojetSel->SetMinEmulMetEt(iCat, 0.);
    monojetSel->SetMinChargedHadronFrac(iCat, 0.2); 
    monojetSel->SetMaxNeutralHadronFrac(iCat, 0.7);
    monojetSel->SetMaxNeutralEmFrac(iCat, 0.7);
  }

  
  //------------------------------------------------------------------------------------------------
  // making the analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel       ->Add(generatorMod);
  generatorMod     ->Add(goodPvMod);
  goodPvMod        ->Add(hltModP);
  // photon regression
  hltModP          ->Add(photonReg);
  
  // simple object id modules
  photonReg        ->Add(SepPUMod); 
  SepPUMod         ->Add(muonId);
  
  muonId           ->Add(eleIdMod);
  eleIdMod	   ->Add(electronCleaning);
  electronCleaning ->Add(merger);
  merger           ->Add(photonIDMod);
  photonIDMod	   ->Add(photonCleaningMod);
  photonCleaningMod->Add(pfTauIDMod);
  pfTauIDMod       ->Add(pfTauCleaningMod);
  pfTauCleaningMod ->Add(pubJet);
  pubJet           ->Add(jetCorr);
  jetCorr          ->Add(metCorrT0T1Shift);
  metCorrT0T1Shift ->Add(jetID);
  jetID            ->Add(jetCleaning);
  jetCleaning      ->Add(monojetSel);
 
  //monojetSel->Add(counter);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseCacher(1);
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimDataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset, 1); // 1 to use smartcache
  else
    d = c->FindDataset(book,skimDataset.Data(),fileset, 1);
  ana->AddDataset(d);
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // skim output
  //------------------------------------------------------------------------------------------------
  TString outputName = TString(outputLabel);
  outputName += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    outputName += TString("_") + TString(fileset);
  
  OutputMod *skimOutput = new OutputMod;
  skimOutput->Drop("*");
  skimOutput->Keep("HLT*");
  
  skimOutput->Keep("MC*");
  skimOutput->Keep("PileupInfo");
  skimOutput->Keep("Rho");
  skimOutput->Keep("EvtSelData");
  skimOutput->Keep("BeamSpot");
  skimOutput->Keep("PrimaryVertexes");
  skimOutput->Keep("PFMet");
  skimOutput->Keep("AKt5PFJets");
  skimOutput->Keep("Electrons");
  skimOutput->Keep("Muons");
  skimOutput->Keep("HPSTaus");
  skimOutput->Keep("Photons");
  // TODO find some object that is named and can hold bit mask info
  skimOutput->AddNewBranch(monojetSel->GetCategoriesName());

  skimOutput->SetMaxFileSize(10 * 1024); // 10 GB - should never exceed
  skimOutput->SetFileName(outputName);
  skimOutput->SetPathName(".");


  //hltModP->Add(monojetSel);
  //counter->Add(skimOutput);
  monojetSel->Add(skimOutput);




  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n",jsonFile.Data());
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(false);

  delete ana; // all modules deleted recursively

  // rename the output file so that condor can see it
  gSystem->Rename("./" + outputName + "_000.root", "./" + outputName + ".root");

  return;
}

