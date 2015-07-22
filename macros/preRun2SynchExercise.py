# Pre-run2 monojet synchronization exercise.
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/DMS13TeVSynch
#
# The objective was to synchronize Bambu produced from MiniAOD with
# other groups who are using MiniAOD directly to produce their ntuples.
# Through the exercise, a problem was found in the Jet Id. The PU
# Jet Id MVA value obtained from MiniAOD jet.userFloat() is a quantity
# calculated from the full PFJet during the production of MiniAOD, and
# is not reproducible from information in the MiniAOD (and thus in Bambu
# files produced from it). Therefore, for full synchronization, a hack
# to read this precalculated MVA values from MiniAOD and store them as
# an array in the Bambu file was introduced. This version of the macro
# only contains on-the-fly JetId calculation (and therefore does not
# agree perfectly with other groups).

from MitAna.TreeMod.bambu import mithep, analysis
import os

mitdata = os.environ['MIT_DATA']

from MitPhysics.Mods.GoodPVFilterMod import goodPVFilterMod
from MitPhysics.Mods.SeparatePileUpMod import separatePileUpMod
from MitPhysics.Mods.ElectronIdMod import electronIdMod
from MitPhysics.Mods.MuonIdMod import muonIdMod
from MitPhysics.Mods.PFTauIdMod import pfTauIdMod
from MitPhysics.Mods.PhotonIdMod import photonIdMod
from MitPhysics.Mods.JetCorrectionMod import jetCorrectionMod
from MitPhysics.Mods.JetIdMod import jetIdMod
from MitPhysics.Mods.MetCorrectionMod import metCorrectionMod

preRun2Sync = mithep.PreRun2SynchExercise(
    VerticesName = mithep.Names.gkPVBrn,
    MetName = metCorrectionMod.GetOutputName(),
    JetsName = jetIdMod.GetOutputName(),
    ElectronsName = electronIdMod.GetOutputName(),
    MuonsName = muonIdMod.GetOutputName(),
    TausName = pfTauIdMod.GetOutputName(),
    PhotonsName = photonIdMod.GetOutputName()
)

analysis.setSequence(
    goodPVFilterMod *
    separatePileUpMod *
    electronIdMod *
    muonIdMod *
    pfTauIdMod *
    photonIdMod *
    jetCorrectionMod *
    jetIdMod *
    metCorrectionMod *
    preRun2Sync
)
