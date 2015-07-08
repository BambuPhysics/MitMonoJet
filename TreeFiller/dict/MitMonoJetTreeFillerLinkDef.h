#ifndef MITMONOJET_TREEFILLER_LINKDEF_H
#define MITMONOJET_TREEFILLER_LINKDEF_H
#include "MitMonoJet/TreeFiller/interface/FillerXlJets.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlMet.h"
#include "MitMonoJet/TreeFiller/interface/FillerXsIsoParticles.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::FillerXlJets+;
#pragma link C++ class mithep::FillerXlMet+;
#pragma link C++ class mithep::FillerXsIsoParticles+;
#endif
