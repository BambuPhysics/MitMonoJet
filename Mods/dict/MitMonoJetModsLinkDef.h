#ifndef MITMONOJET_MODS_LINKDEF_H
#define MITMONOJET_MODS_LINKDEF_H
#include "MitMonoJet/Mods/interface/HltEvtSelMod.h"
#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"
#include "MitMonoJet/Mods/interface/BoostedVTreeWriter.h"
#include "MitMonoJet/Mods/interface/DMSTreeWriter.h"
#include "MitMonoJet/Mods/interface/SkimJetsMod.h"
#include "MitMonoJet/Mods/interface/MetAnalysisMod.h"
#include "MitMonoJet/Mods/interface/PreRun2SynchExercise.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::HltEvtSelMod+;
#pragma link C++ class mithep::MonoJetEventHlt+;
#pragma link C++ class mithep::MonoJetTreeWriter+;
#pragma link C++ class mithep::BoostedVTreeWriter+;
#pragma link C++ class mithep::DMSTreeWriter+;
#pragma link C++ class mithep::SkimJetsMod+;
#pragma link C++ class mithep::MetAnalysisMod+;
#pragma link C++ class mithep::PreRun2SynchExercise+;
#endif
