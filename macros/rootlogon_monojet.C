{
  // ------------------------------------------------------------------------------------------------
  //
  // Root logon for MIT analysis
  // 
  // It is important to note that we are not loading any shared library explicitely but we rely on
  // the compiler for the macro you execute to do that for us. In this way the minimum amount of
  // libraries are being loaded and the resulting shared library file has all libraries created such
  // that they are not hardmounted to a directory. This is un-important if you are working on a
  // cluster with common home (ex. NFS) but if you want to transplant the job into another system it
  // becomes crucial.
  //
  // ------------------------------------------------------------------------------------------------
  // To see what is happening behind the scenes use this
  gDebug = kFALSE;

  // Maintain the direcories and libraries we need to use for root internal compilation
  TString str = gSystem->GetMakeSharedLib();
  str += TString(" -L$CMSSW_BASE/lib/slc6_amd64_gcc491");
  str += TString(" -L$CMSSW_BASE/external/slc6_amd64_gcc491/lib");
  str += TString(" -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/fastjet/3.1.0-odfocd/lib");
  str += TString(" -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/fastjet-contrib/1.014-odfocd/lib");
  str += " -lMitCommonDataFormats";
  str += TString(" -lMitAnaCatalog -lMitAnaDataCont -lMitAnaDataTree -lMitAnaDataUtil");
  str += TString(" -lMitAnaPhysicsMod -lMitAnaTAM -lMitAnaTreeMod -lMitAnaUtils");
  str += TString(" -lMitAnaValidation");
  str += TString(" -lMitPhysicsMods -lMitPhysicsSelMods -lMitPhysicsUtils");
  str += TString(" -lMitMonoJetMods -lMitMonoJetSelMods -lMitMonoJetTreeFiller");
  str += TString(" -lfastjet -lfastjettools -lfastjetcontribfragile -lqjets");
  str += TString(" -lMitPlotsStyle -lMitPlotsInput -lMitPlotsPlot");
  gSystem->SetMakeSharedLib(str);

  if (gDebug)
    printf(" Shared libraries: %s\n",gSystem->GetMakeSharedLib());

  // Make sure we have all include files
  gSystem->AddIncludePath("-I$EXTERNAL/include");
  gSystem->AddIncludePath("-I$EXTERNAL/Qjets");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/MitAna/macros");
  
  gInterpreter->AddIncludePath(TString("/cvmfs/cms.cern.ch/")+TString(gSystem->Getenv("SCRAM_ARCH"))+
			       TString("/lcg/roofit/5.32.00/include"));
  gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_BASE"))+TString("/src/"));
  gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_RELEASE_BASE"))+TString("/src/"));
  gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_BASE"))+
			       TString("/src/MitAna/TreeMod/interface"));
  
  gROOT->SetMacroPath(TString(gROOT->GetMacroPath())+":"+
		      TString("$CMSSW_BASE/src/MitPlots/macros"));
}
