#!/bin/bash

$CMSSW_BASE/src/MitMonoJet/bin/installQjets.sh

if [ -e $HOME/cms/root/rootlogon.C ]
then
  cp $HOME/cms/root/rootlogon.C $HOME/cms/root/rootlogon.C.last
  echo " File"
  echo "$HOME/cms/root/rootlogon.C"
  echo " copied to"
  echo "$HOME/cms/root/rootlogon.C.last"
fi

cp $CMSSW_BASE/src/MitMonoJet/macros/rootlogon_monojet.C $HOME/cms/root/rootlogon.C

$CMSSW_BASE/src/MitCommon/bin/genDict.sh MitMonoJet/{Mods,SelMods,TreeFiller,Utils}