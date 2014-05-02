//--------------------------------------------------------------------------------------------------
// $Id: FillerXlJets.h,v 1.9 2011/03/01 17:27:22 mzanetti Exp $
//
// FillerXlJets
//
// This module process a collection of input jet, compute the substrucure
// and fill two output collections of fXlFatJets and relative fXlSubJets
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_TREEFILLER_FILLERXLJETS_H
#define MITMONOJET_TREEFILLER_FILLERXLJETS_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "MitMonoJet/DataTree/interface/XlFatJetFwd.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"
#include "MitMonoJet/DataTree/interface/XlSubJetFwd.h"
#include "MitMonoJet/DataTree/interface/XlSubJet.h"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

namespace mithep
{
  class FillerXlJets : public BaseMod
  {
    public:
      FillerXlJets(const char *name = "FillerXlJets",
                   const char *title = "XlJets Filler module");
      ~FillerXlJets();

      void IsData(Bool_t b)                { fIsData = b;         }
      void FillVSubJets(Bool_t b)          { fFillVSubJets = b;   }
      void FillTopSubJets(Bool_t b)        { fFillTopSubJets = b; }
      void SetBtaggingOn(Bool_t b)         { fBTaggingActive = b; }
      void SetfQGTaggingOn(Bool_t b)       { fQGTaggingActive = b;}
      void PublishOutput(Bool_t b)         { fPublishOutput = b;  }

      void SetProcessNJets(UInt_t n)       { fProcessNJets = n;  } 
      
      void SetJetsName(const char *n)      { fJetsName = n;       }
      void SetJetsFromBranch(Bool_t b)     { fJetsFromBranch = b; }

      void SetFatJetsName(const char *n)   { fXlFatJetsName = n;  }
      void SetSubJetsName(const char *n)   { fXlSubJetsName = n;  }
       
      void SetPruningOn(Bool_t b)          { fPrune = b;          }
      void SetFilteringOn(Bool_t b)        { fFilter = b;         }
      void SetTrimmingOn(Bool_t b)         { fTrim = b;           }
      
      void SetPruneZCut(double d)          { fPruneZCut = d;      }
      void SetPruneDistCut(double d)       { fPruneDistCut = d;   }
      void SetFilterN(int n)               { fFilterN = n;        }
      void SetFilterRad(double d)          { fFilterRad = d;      }
      void SetTrimRad(double d)            { fTrimRad = d;        }
      void SetTrimPtFrac(double d)         { fTrimPtFrac = d;     }
      void SetConeSize(double d)           { fConeSize = d;       }
 
    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
 
      void FillfXlFatJets(std::vector<fastjet::PseudoJet> &fjFatJets);
      void FillfXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets, XlFatJet *pFatJet,
                         XlSubJet::ESubJetType t);
 
    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fFillVSubJets;                //=true if V-subjets are stored (2-prom structure)
      Bool_t fFillTopSubJets;              //=true if top-subjets are stored (3-prom structure)
      Bool_t fBTaggingActive;              //=true if BTagging info is filled
      Bool_t fQGTaggingActive;             //=true if QGTagging info is filled
      Bool_t fPublishOutput;               //=true if output collection are published

      UInt_t  fProcessNJets;               //number of input jets processed by fastjet

      TString fJetsName;                   //(i) name of input jets
      Bool_t fJetsFromBranch;              //are input jets from Branch?
      const JetCol *fJets;                 //input jets

      TString fPfCandidatesName;           //(i) name of PF candidates coll
      Bool_t fPfCandidatesFromBranch;      //are PF candidates from Branch?
      const PFCandidateCol *fPfCandidates; //particle flow candidates coll handle
 
      TString fXlFatJetsName;              //name of output fXlFatJets collection
      XlFatJetArr *fXlFatJets;             //array of fXlFatJets
      TString fXlSubJetsName;              //name of output fXlSubJets collection
      XlSubJetArr *fXlSubJets;             //array of fXlSubJets
      
      // Objects from fastjet we want to use
      Bool_t fPrune;                       //apply pruning?
      Bool_t fFilter;                      //apply filtering?
      Bool_t fTrim;                        //apply trimming?
      fastjet::Pruner *fPruner;
      fastjet::Filter *fFilterer;
      fastjet::Filter *fTrimmer;           //no this is not a typo trimmers belong to fastjet Filter class
      double fPruneZCut;                   //pruning Z cut
      double fPruneDistCut;                //pruning distance cut 
      int fFilterN;                        //number of subjets after filtering
      double fFilterRad;                   //filtered subjets radius
      double fTrimRad;                     //trimmed subjet radius
      double fTrimPtFrac;                  //trimmed subjet pt fraction
      double fConeSize;                    //fastjet clustering radius
      fastjet::JetDefinition *fCAJetDef;   //fastjet clustering definition
      fastjet::GhostedAreaSpec *fActiveArea;
      fastjet::AreaDefinition *fAreaDefinition;

      ClassDef(FillerXlJets, 0)            //XlJets, Fat and Sub, filler      
  };
}
#endif
