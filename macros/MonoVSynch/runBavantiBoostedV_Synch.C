#include <TSystem.h>
#include <TProfile.h>
#include "MitCommon/Utils/interface/Utils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/SkimMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitMonoJet/SelMods/interface/BoostedVAnalysisMod.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlJets.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlFatJets.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlMet.h"
#include "MitMonoJet/TreeFiller/interface/FillerXsIsoParticles.h"
#include "MitMonoJet/Mods/interface/FastJetMod.h"
#include "MitMonoJet/Mods/interface/SkimJetsMod.h"

TString getCatalogDir(const char* dir);
TString getJsonFile(const char* dir);

//--------------------------------------------------------------------------------------------------
void runBavantiBoostedV_Synch
                       (const char *fileset    = "0000",
                        const char *skim       = "noskim",
                        //const char *dataset    = "s12-ttj-sl-v1-v7c",     
                        const char *dataset    = "GJets_HT-400ToInf_8TeV-madgraph+Summer12_DR53X-PU_S10_START53_V7A-v1+AODSIM",                        
                        const char *book       = "t2mit/filefi/032",
                        const char *catalogDir = "/home/cmsprod/catalog",
                        const char *outputName = "boostedv",
                        int         nEvents    = 10)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  TString cataDir  = getCatalogDir(catalogDir);
  TString mitData  = Utils::GetEnv("MIT_DATA");
  TString json     = Utils::GetEnv("MIT_PROD_JSON");
  TString jsonFile = getJsonFile("/home/cmsprod/cms/json");
  Bool_t  isData   = (json.CompareTo("~") != 0);
  // only for synchrionization purposes
  isData = kFALSE;
  printf("\n Initialization worked. Data?: %d\n\n",isData);

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  // debugging config
  using namespace mithep;
  gDebugMask  = (Debug::EDebugMask) (Debug::kGeneral | Debug::kTreeIO);
  gDebugLevel = 3;


  // Caching and how
  Int_t local = 1, cacher = 1;

  // local =   0 - as is,
  //           1 - /mt/hadoop  (MIT:SmartCache - preload one-by-one)
  //           2 - /mnt/hadoop (MIT:SmartCache - preload complete fileset)
  //           3 - ./          (xrdcp          - preload one-by-one)
  // cacher =  0 - no file by file caching
  //           1 - file by file caching on

  //------------------------------------------------------------------------------------------------
  // set up information for master module
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted
  
  // only select on run- and lumisection numbers when valid json file present
  if (json.CompareTo("~") != 0 && json.CompareTo("-") != 0) {
    printf(" runBoostedV() - adding jsonFile: %s\n",jsonFile.Data());
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if (json.CompareTo("-") == 0) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }
  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseCacher(cacher);
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(cataDir.Data());
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset,local);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset,local);
  ana->AddDataset(d);
  
  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  TString ntupleFile = rootFile + TString("_ntuple");
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModP = new HLTMod("HLTModP");
  hltModP->SetBitsName("HLTBits");
  hltModP->SetTrigObjsName("HltObjsMonoJet");
  hltModP->SetAbortIfNotAccepted(isData);
  hltModP->SetPrintTable(kFALSE);

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
                                         "HLT_MET120_HBHENoiseCleaned_v2"};

  for (int i=0; i<nMjtTrigs; i++)
    hltModP->AddTrigger(TString("!+"+monoJetTriggers[i]),0,999999);

  // top semileptonic triggers
  const int nTopTrigs = 14;
  TString topTriggers[nTopTrigs] = { "HLT_IsoMu15_v2",
                                     "HLT_IsoMu24_v2",
                                     "HLT_IsoMu17_v6",
                                     "HLT_IsoMu17_v8",
                                     "HLT_IsoMu17_v9",
                                     "HLT_IsoMu17_eta2p1_v1",
                                     "HLT_IsoMu24_v8", 
                                     "HLT_IsoMu24_eta2p1_v3", 
                                     "HLT_IsoMu24_eta2p1_v6", 
                                     "HLT_IsoMu24_eta2p1_v7", 
                                     "HLT_IsoMu24_eta2p1_v12", 
                                     "HLT_IsoMu24_eta2p1_v13", 
                                     "HLT_IsoMu24_eta2p1_v14", 
                                     "HLT_IsoMu24_eta2p1_v15"};

  for (int i=0; i<nTopTrigs; i++)
    hltModP->AddTrigger(TString("!+"+topTriggers[i]),0,999999);

  // photon triggers
  const int nPhotonTrigs = 6;
  TString photonTriggers[nPhotonTrigs] = { "HLT_Photon135_v4",
                                           "HLT_Photon135_v5",
                                           "HLT_Photon135_v6",
                                           "HLT_Photon135_v7",
                                           "HLT_Photon150_v3",
                                           "HLT_Photon150_v4"};

  for (int i=0; i<nPhotonTrigs; i++)
    hltModP->AddTrigger(TString("!+"+photonTriggers[i]),0,999999);

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

  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* sepPuMod = new SeparatePileUpMod;
  //sepPuMod->SetUseAllVerteces(kFALSE);
  //sepPuMod->SetVertexName("OutVtxCiC");
  sepPuMod->SetPFNoPileUpName("pfnopileupcands");
  sepPuMod->SetPFPileUpName("pfpileupcands");
  sepPuMod->SetCheckClosestZVertex(kFALSE);

  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4.0);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  goodPVFilterMod->SetIsMC(!isData);
  goodPVFilterMod->SetVertexesName("PrimaryVertexes");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------

  //-----------------------------------
  // Lepton Selection
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod->SetPtMin(10.);
  eleIdMod->SetEtaMax(2.5);
  eleIdMod->SetApplyEcalFiducial(kTRUE);
  eleIdMod->SetIDType("VBTFWorkingPoint95Id");
  eleIdMod->SetIsoType("PFIso");
  eleIdMod->SetApplyConversionFilterType1(kTRUE);
  eleIdMod->SetApplyConversionFilterType2(kFALSE);
  eleIdMod->SetChargeFilter(kFALSE);
  eleIdMod->SetApplyD0Cut(kTRUE);
  eleIdMod->SetApplyDZCut(kTRUE);
  eleIdMod->SetWhichVertex(0);
  eleIdMod->SetNExpectedHitsInnerCut(0);
  eleIdMod->SetGoodElectronsName("GoodElectronsBS");
  eleIdMod->SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);

  MuonIDMod *muonIdIsoMod = new MuonIDMod;
  muonIdIsoMod->SetOutputName("HWWMuons");
  muonIdIsoMod->SetIntRadius(0.0);
  muonIdIsoMod->SetClassType("GlobalTracker");
  muonIdIsoMod->SetIDType("WWMuIdV4");
  muonIdIsoMod->SetIsoType("IsoRingsV0_BDTG_Iso");
  muonIdIsoMod->SetApplyD0Cut(kTRUE);
  muonIdIsoMod->SetD0Cut(0.02);
  muonIdIsoMod->SetApplyDZCut(kTRUE);
  muonIdIsoMod->SetDZCut(0.1);
  muonIdIsoMod->SetWhichVertex(0);
  muonIdIsoMod->SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  muonIdIsoMod->SetPtMin(10.);
  muonIdIsoMod->SetEtaCut(2.4);

  MuonIDMod *muonIdLooseMod = new MuonIDMod;
  muonIdLooseMod->SetOutputName("POGMuons");
  muonIdLooseMod->SetClassType("GlobalTracker");
  muonIdLooseMod->SetIDType("NoId");
  muonIdLooseMod->SetIsoType("NoIso");
  muonIdLooseMod->SetApplyD0Cut(true);
  muonIdLooseMod->SetD0Cut(0.2);
  muonIdLooseMod->SetApplyDZCut(true);
  muonIdLooseMod->SetDZCut(0.5);
  muonIdLooseMod->SetPtMin(10.);
  muonIdLooseMod->SetEtaCut(2.4);
  
  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  electronCleaning->SetCleanMuonsName(muonIdLooseMod->GetOutputName());
  electronCleaning->SetGoodElectronsName(eleIdMod->GetOutputName());
  electronCleaning->SetCleanElectronsName("CleanElectrons");

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName(muonIdLooseMod->GetOutputName());
  merger->SetElectronsName(electronCleaning->GetOutputName());
  merger->SetMergedName("MergedLeptons");

  //-----------------------------------
  // Photon Egamma Medium ID
  //-----------------------------------
  PhotonIDMod *photonIdMod = new PhotonIDMod;
  photonIdMod->SetPtMin(10.0);
  photonIdMod->SetOutputName("GoodPhotons");
  photonIdMod->SetIDType("EgammaMedium");
  photonIdMod->SetIsoType("NoIso");
  photonIdMod->SetApplyElectronVeto(kTRUE);
  photonIdMod->SetApplyPixelSeed(kFALSE);
  photonIdMod->SetApplyConversionId(kFALSE);
  photonIdMod->SetApplyFiduciality(kTRUE);
  photonIdMod->SetIsData(isData);
  photonIdMod->SetPhotonsFromBranch(kTRUE);

  PhotonCleaningMod *photonCleaningMod = new PhotonCleaningMod;
  photonCleaningMod->SetCleanElectronsName(electronCleaning->GetOutputName());
  photonCleaningMod->SetGoodPhotonsName(photonIdMod->GetOutputName());
  photonCleaningMod->SetCleanPhotonsName("CleanPhotons");

  PFTauIDMod *pftauIdMod = new PFTauIDMod;
  pftauIdMod->SetPFTausName("HPSTaus");
  pftauIdMod->SetIsLooseId(kTRUE);
  pftauIdMod->SetIsHPSSel(kTRUE); // to get >= 5_3_14 samples running
  pftauIdMod->SetPtMin(15); // to loosen the WP
  
  PFTauCleaningMod *pftauCleaningMod = new PFTauCleaningMod;
  pftauCleaningMod->SetGoodPFTausName(pftauIdMod->GetGoodPFTausName());
  pftauCleaningMod->SetCleanMuonsName(muonIdLooseMod->GetOutputName());
  pftauCleaningMod->SetCleanElectronsName(electronCleaning->GetOutputName());

  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt4PFJets");
  pubJet->SetOutputName("PubAKt4PFJets");
  
  FastJetMod *pubFastJet = new FastJetMod;
  pubFastJet->SetPfCandidatesName("pfnopileupcands");
  pubFastJet->SetPfCandidatesFromBranch(kFALSE);
  pubFastJet->SetConeSize(0.8);

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if (isData){ 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data());
  }                                                                                      
  else {                                                                                 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    

  JetCorrectionMod *fatJetCorr = new JetCorrectionMod;
  if (isData){ 
    fatJetCorr->AddCorrectionFromFile((mitData+TString("/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt")).Data()); 
    fatJetCorr->AddCorrectionFromFile((mitData+TString("/FT53_V21A_AN6_L2Relative_AK7PFchs.txt")).Data()); 
    fatJetCorr->AddCorrectionFromFile((mitData+TString("/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt")).Data()); 
    fatJetCorr->AddCorrectionFromFile((mitData+TString("/FT53_V21A_AN6_L2L3Residual_AK7PFchs.txt")).Data());
  }                                                                                      
  else {                                                                                 
    fatJetCorr->AddCorrectionFromFile((mitData+TString("/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt")).Data()); 
    fatJetCorr->AddCorrectionFromFile((mitData+TString("/FT53_V21A_AN6_L2Relative_AK7PFchs.txt")).Data()); 
    fatJetCorr->AddCorrectionFromFile((mitData+TString("/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt")).Data()); 
  }
  fatJetCorr->SetInputName(pubFastJet->GetOutputJetsName());
  fatJetCorr->SetCorrectedName("CorrectedFatJets");    
        
  JetIDMod *jetId = new JetIDMod;
  jetId->SetInputName(jetCorr->GetOutputName());
  jetId->SetPtCut(30.0);
  jetId->SetEtaMaxCut(2.5);
  jetId->SetJetEEMFractionMinCut(0.00);
  jetId->SetOutputName("GoodJets");
  jetId->SetApplyBetaCut(kFALSE);
  jetId->SetApplyMVACut(kTRUE);

  JetIDMod *fatJetId = new JetIDMod;
  fatJetId->SetInputName(fatJetCorr->GetOutputName());
  fatJetId->SetPtCut(30.0);
  fatJetId->SetEtaMaxCut(2.5);
  fatJetId->SetJetEEMFractionMinCut(0.00);
  fatJetId->SetOutputName("GoodFatJets");
  fatJetId->SetApplyBetaCut(kFALSE);
  fatJetId->SetApplyMVACut(kFALSE);
  fatJetId->SetApplyMVACHS(kFALSE);

  JetCleaningMod *jetCleaning = new JetCleaningMod;
  jetCleaning->SetCleanElectronsName(electronCleaning->GetOutputName());
  jetCleaning->SetCleanMuonsName(muonIdIsoMod->GetOutputName());
  jetCleaning->SetCleanPhotonsName(photonCleaningMod->GetOutputName());
  jetCleaning->SetApplyPhotonRemoval(kTRUE);
  jetCleaning->SetApplyTauRemoval(kFALSE);
  jetCleaning->SetGoodJetsName(jetId->GetOutputName());
  jetCleaning->SetCleanJetsName("CleanJets");

  JetCleaningMod *fatJetCleaning = new JetCleaningMod;
  fatJetCleaning->SetCleanElectronsName(electronCleaning->GetOutputName());
  fatJetCleaning->SetCleanMuonsName(muonIdIsoMod->GetOutputName());
  fatJetCleaning->SetCleanPhotonsName(photonCleaningMod->GetOutputName());
  fatJetCleaning->SetMinDeltaRToElectron(0.5);
  fatJetCleaning->SetMinDeltaRToMuon(0.5);
  fatJetCleaning->SetMinDeltaRToPhoton(0.5);
  fatJetCleaning->SetApplyPhotonRemoval(kTRUE);
  fatJetCleaning->SetApplyTauRemoval(kFALSE);  
  fatJetCleaning->SetGoodJetsName(fatJetId->GetOutputName());
  fatJetCleaning->SetCleanJetsName("CleanFatJets");

  //------------------------------------------------------------------------------------------------
  // select events with a monojet topology
  //------------------------------------------------------------------------------------------------
  BoostedVAnalysisMod *jetplusmet = new BoostedVAnalysisMod("MonoJetSelector");
  jetplusmet->SetFatJetsName(fatJetCleaning->GetOutputName()); //identified fat jets
  jetplusmet->SetFatJetsFromBranch(kFALSE);
  jetplusmet->SetJetsName(jetCleaning->GetOutputName()); //identified jets
  jetplusmet->SetJetsFromBranch(kFALSE);
  jetplusmet->SetElectronsName(electronCleaning->GetOutputName());
  jetplusmet->SetElectronsFromBranch(kFALSE);
  jetplusmet->SetMuonsName(muonIdLooseMod->GetOutputName());
  jetplusmet->SetMuonsFromBranch(kFALSE);
  jetplusmet->SetPhotonsName(photonCleaningMod->GetOutputName());
  jetplusmet->SetPhotonsFromBranch(kFALSE);
  jetplusmet->ApplyResolvedPresel(kTRUE); 
  jetplusmet->ApplyTopPresel(kFALSE); 
  jetplusmet->ApplyWlepPresel(kTRUE);
  jetplusmet->ApplyZeePresel(kTRUE);
  jetplusmet->ApplyZmmPresel(kTRUE);
  jetplusmet->ApplyMetPresel(kTRUE);
  jetplusmet->ApplyVbfPresel(kFALSE);
  jetplusmet->ApplyGjetPresel(kTRUE);
  jetplusmet->SetMinFatJetPt(200);
  jetplusmet->SetMinTagJetPt(110);
  jetplusmet->SetMinMet(150);    
  jetplusmet->SetMinPhotonPt(160);    
  jetplusmet->SkipEvents(kFALSE);    

  //------------------------------------------------------------------------------------------------
  // prepare the extended MVA met 
  //------------------------------------------------------------------------------------------------
  FillerXlMet *extendedMetFiller = new FillerXlMet();
  extendedMetFiller->SetIsData(isData);
  extendedMetFiller->SetJetsFromBranch(kFALSE);
  extendedMetFiller->SetJetsName(pubJet->GetOutputName());
  extendedMetFiller->SetMuonsFromBranch(kTRUE);
  //extendedMetFiller->SetMuonsName(muonIdLooseMod->GetOutputName());
  extendedMetFiller->SetElectronsFromBranch(kTRUE);
  //extendedMetFiller->SetElectronsName(electronCleaning->GetOutputName());
  extendedMetFiller->SetTausFromBranch(kTRUE);
  extendedMetFiller->SetTausName("HPSTaus");
  extendedMetFiller->SetPhotonsFromBranch(kTRUE);
  //extendedMetFiller->SetPhotonsName(photonCleaningMod->GetOutputName());
  extendedMetFiller->SetPVFromBranch(kFALSE);
  extendedMetFiller->SetPVName(goodPVFilterMod->GetOutputName());
  extendedMetFiller->SetXlMetName("PFMetMVA");     

  //------------------------------------------------------------------------------------------------
  // prepare the extended jets with QG/color pull information
  //------------------------------------------------------------------------------------------------
  FillerXlJets *boostedJetsFiller = new FillerXlJets;  
  boostedJetsFiller->SetJetsName(jetCleaning->GetOutputName());
  boostedJetsFiller->SetJetsFromBranch(kFALSE);
  
  //------------------------------------------------------------------------------------------------
  // prepare the extended jets with substructure information
  //------------------------------------------------------------------------------------------------
  FillerXlFatJets *boostedFatJetsFiller = new FillerXlFatJets;  
  boostedFatJetsFiller->FillTopSubJets(kFALSE);
  boostedFatJetsFiller->SetJetsName(fatJetCleaning->GetOutputName());
  boostedFatJetsFiller->SetJetsFromBranch(kFALSE);
  boostedFatJetsFiller->SetConeSize(0.8);      

  //------------------------------------------------------------------------------------------------
  // prepare the reduced isolated particles
  //------------------------------------------------------------------------------------------------
  FillerXsIsoParticles *boostedXsIsoParticlesFiller = new FillerXsIsoParticles;  
  boostedXsIsoParticlesFiller->FillXsMuons(kTRUE);
  boostedXsIsoParticlesFiller->FillXsElectrons(kTRUE);
  boostedXsIsoParticlesFiller->FillXsTaus(kTRUE);
  boostedXsIsoParticlesFiller->FillXsPhotons(kTRUE);
  boostedXsIsoParticlesFiller->SetMuonsName(muonIdLooseMod->GetOutputName());
  boostedXsIsoParticlesFiller->SetMuonsFromBranch(kFALSE);
  boostedXsIsoParticlesFiller->SetIsoMuonsName(muonIdIsoMod->GetOutputName());
  boostedXsIsoParticlesFiller->SetIsoMuonsFromBranch(kFALSE);
  boostedXsIsoParticlesFiller->SetElectronsName(eleIdMod->GetOutputName());
  boostedXsIsoParticlesFiller->SetElectronsFromBranch(kFALSE);
  boostedXsIsoParticlesFiller->SetTausName(pftauCleaningMod->GetOutputName());
  boostedXsIsoParticlesFiller->SetTausFromBranch(kFALSE);
  boostedXsIsoParticlesFiller->SetPhotonsName(photonCleaningMod->GetOutputName());
  boostedXsIsoParticlesFiller->SetPhotonsFromBranch(kFALSE);

  //------------------------------------------------------------------------------------------------
  // keep the skimmed collections for further usage
  //------------------------------------------------------------------------------------------------
  //SkimMod<PFCandidate> *skmPFCandidates = new SkimMod<PFCandidate>;
  //skmPFCandidates->SetBranchName(Names::gkPFCandidatesBrn);
  //skmPFCandidates->SetPublishArray(kTRUE);

  SkimJetsMod *skmJets = new SkimJetsMod;
  skmJets->SetBranchName(jetCleaning->GetOutputName());
  skmJets->SetColFromBranch(kFALSE);
  skmJets->SetColMarkFilter(kFALSE);
  skmJets->SetPublishArray(kTRUE);

  SkimMod<TriggerObject> *skmTrigger = new SkimMod<TriggerObject>;
  skmTrigger->SetBranchName(hltModP->GetOutputName());
  skmTrigger->SetColFromBranch(kFALSE);
  skmTrigger->SetColMarkFilter(kFALSE);
  skmTrigger->SetPublishArray(kTRUE);
    
  //------------------------------------------------------------------------------------------------
  // save all this in an output ntuple
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->SetUseBrDep(kFALSE);
  outMod->SetKeepTamBr(kFALSE);
  outMod->SetFileName(ntupleFile);
  outMod->SetMaxFileSize(4096);
  outMod->Drop("*");
  outMod->Keep(Names::gkMCEvtInfoBrn);
  outMod->Keep(Names::gkMCPartBrn);
  outMod->Keep(Names::gkPVBeamSpotBrn);
  outMod->Keep(Names::gkPileupInfoBrn);
  outMod->Keep(Names::gkPileupEnergyDensityBrn);
  outMod->Keep("PFMet");
  outMod->AddNewBranch("XlEvtSelData");
  //outMod->AddNewBranch(TString("Skm") + Names::gkPFCandidatesBrn);
  outMod->AddNewBranch(TString("Skm") + hltModP->GetOutputName());
  outMod->AddNewBranch("PFMetMVA");
  outMod->AddNewBranch("XlJets");
  outMod->AddNewBranch("XlFatJets");
  outMod->AddNewBranch("XlSubJets");
  outMod->AddNewBranch("XsMuons");
  outMod->AddNewBranch("XsElectrons");
  outMod->AddNewBranch("XsTaus");
  outMod->AddNewBranch("XsPhotons");
  
  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  runLumiSel                 ->Add(goodPVFilterMod);
  goodPVFilterMod            ->Add(hltModP);
  hltModP                    ->Add(sepPuMod);
  sepPuMod                   ->Add(muonIdIsoMod);
  muonIdIsoMod               ->Add(muonIdLooseMod);
  muonIdLooseMod             ->Add(eleIdMod);
  eleIdMod                   ->Add(electronCleaning);
  electronCleaning           ->Add(merger);
  merger                     ->Add(photonIdMod);
  photonIdMod                ->Add(photonCleaningMod);
  photonCleaningMod          ->Add(pftauIdMod);
  pftauIdMod                 ->Add(pftauCleaningMod);
  pftauCleaningMod           ->Add(pubJet);
  pubJet                     ->Add(jetCorr);
  jetCorr                    ->Add(jetId);
  jetId                      ->Add(jetCleaning);
  jetCleaning                ->Add(pubFastJet);
  pubFastJet                 ->Add(fatJetCorr);
  fatJetCorr                 ->Add(fatJetId);
  fatJetId                   ->Add(fatJetCleaning);
  fatJetCleaning             ->Add(jetplusmet);
  jetplusmet                 ->Add(extendedMetFiller);
  extendedMetFiller          ->Add(boostedJetsFiller);
  boostedJetsFiller          ->Add(boostedFatJetsFiller);
  boostedFatJetsFiller       ->Add(boostedXsIsoParticlesFiller);
  boostedXsIsoParticlesFiller->Add(skmJets);
  skmJets                    ->Add(skmTrigger);
  skmTrigger                 ->Add(outMod);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n",jsonFile.Data());
  printf("\n Rely on Catalog: %s\n",cataDir.Data());
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output:   %s\n",rootFile.Data());  
  printf("\n Ntuple output: %s\n\n",(ntupleFile + TString(".root")).Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}

//--------------------------------------------------------------------------------------------------
TString getCatalogDir(const char* dir)
{
  TString cataDir = TString("./catalog");
  Long_t *id=0,*size=0,*flags=0,*mt=0;

  printf(" Try local catalog first: %s\n",cataDir.Data());
  if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
    cataDir = TString(dir);
    if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
      printf(" Requested local (./catalog) and specified catalog do not exist. EXIT!\n");
      return TString("");
    }
  }
  else {
    printf(" Local catalog exists: %s using this one.\n",cataDir.Data()); 
  }

  return cataDir;
}

//--------------------------------------------------------------------------------------------------
TString getJsonFile(const char* dir)
{
  TString jsonDir  = TString("./json");
  TString json     = Utils::GetEnv("MIT_PROD_JSON");
  Long_t *id=0,*size=0,*flags=0,*mt=0;

  printf(" Try local json first: %s\n",jsonDir.Data());
  if (gSystem->GetPathInfo(jsonDir.Data(),id,size,flags,mt) != 0) {
    jsonDir = TString(dir);
    if (gSystem->GetPathInfo(jsonDir.Data(),id,size,flags,mt) != 0) {
      printf(" Requested local (./json) and specified json directory do not exist. EXIT!\n");
      return TString("");
    }
  }
  else {
    printf(" Local json directory exists: %s using this one.\n",jsonDir.Data()); 
  }

  // Construct the full file name
  TString jsonFile = jsonDir + TString("/") + json;
  if (gSystem->GetPathInfo(jsonFile.Data(),id,size,flags,mt) != 0) {
    printf(" Requested jsonfile (%s) does not exist. EXIT!\n",jsonFile.Data());
    return TString("");
  }
  else {
    printf(" Requested jsonfile (%s) exist. Moving on now!\n",jsonFile.Data());
  }

  return jsonFile;
}
