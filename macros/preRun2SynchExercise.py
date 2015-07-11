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

goodPV = mithep.GoodPVFilterMod(
    MinVertexNTracks = 0,
    MinNDof = 4,
    MaxAbsZ = 24.,
    MaxRho = 2.,
    IsMC = True,
    VertexesName = mithep.Names.gkPVBrn
)

pfPU = mithep.SeparatePileUpMod(
    PFNoPileUpName = "pfNoPU",
    PFPileUpName = "pfPU",
    CheckClosestZVertex = False
)

eleId = mithep.ElectronIdMod(
    OutputName = 'VetoElectrons',
    IdType = mithep.ElectronTools.kPhys14Veto,
    IsoType = mithep.ElectronTools.kPhys14VetoIso,
    ApplyEcalFiducial = True,
    WhichVertex = 0,
    PtMin = 10.,
    EtaMax = 2.5,
    ConversionsName = 'Conversions'
)

muId = mithep.MuonIdMod(
    OutputName = 'LooseMuons',
    IdType = mithep.MuonTools.kNoId,
    IsoType = mithep.MuonTools.kPFIsoBetaPUCorrected,
    PFNoPileupCandidatesName = 'pfNoPU',
    PFPileupCandidatesName = 'pfPU',
    PtMin = 10.,
    EtaMax = 2.4
)

tauId = mithep.PFTauIdMod(
    OutputName = 'LooseTaus',
    PtMin = 18.,
    EtaMax = 2.3
)
tauId.AddDiscriminator(mithep.PFTau.kDiscriminationByDecayModeFindingNewDMs)

photonId = mithep.PhotonIdMod(
    OutputName = 'LoosePhotons',
    IdType = mithep.PhotonTools.kPhys14Loose,
    IsoType = mithep.PhotonTools.kPhys14Loose,
    PtMin = 15.,
    EtaMax = 2.5
)

jetCorr = mithep.JetCorrectionMod(
    InputName = 'AKt4PFJetsCHS',
    CorrectedJetsName = 'AKt4PFJetsCHSL1L2L3',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll
)
jetCorr.AddCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L1FastJet_AK4PFchs.txt")
jetCorr.AddCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L2Relative_AK4PFchs.txt")
jetCorr.AddCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L3Absolute_AK4PFchs.txt")

# UseClassicBetaForMVA: set to true for MiniAOD input. Makes the MVA value agree
# with what CMSSW would give if Jet Id was calculated on MiniAOD. This is not the
# same as the pre-calculated MVA value stored in MiniAOD and accessible with
# jet.userFloat().
jetId = mithep.JetIdMod(
    InputName = jetCorr.GetOutputName(),
    OutputName = 'CleanedJets',
    ApplyPFLooseId = True,
    MVATrainingSet = mithep.JetIDMVA.k53BDTCHSFullPlusRMS,
    MVACutWP = mithep.JetIDMVA.kLoose,
    MVAWeightsFile = mitdata + '/TMVAClassification_5x_BDT_chsFullPlusRMS.weights.xml',
    MVACutsFile = mitdata + '/jetIDCuts_121221.dat',
    UseClassicBetaForMVA = True,
    PtMin = 30.,
    EtaMax = 2.5
)

metCorr = mithep.MetCorrectionMod(
    InputName = 'PFMet',
    OutputName = 'PFType1CorrectedMet',
    JetsName = 'AKt4PFJets',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll,
    MaxEMFraction = 0.9,
    SkipMuons = True
)
metCorr.ApplyType0(False)
metCorr.ApplyType1(True)
metCorr.ApplyShift(False)
metCorr.AddJetCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L1FastJet_AK4PF.txt")
metCorr.AddJetCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L2Relative_AK4PF.txt")
metCorr.AddJetCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L3Absolute_AK4PF.txt")
metCorr.IsData(False)

preRun2Sync = mithep.PreRun2SynchExercise(
    VerticesName = mithep.Names.gkPVBrn,
    MetName = metCorr.GetOutputName(),
    JetsName = jetId.GetOutputName(),
    ElectronsName = eleId.GetOutputName(),
    MuonsName = muId.GetOutputName(),
    TausName = tauId.GetOutputName(),
    PhotonsName = photonId.GetOutputName()
)

analysis.setSequence(
    goodPV *
    pfPU *
    eleId *
    muId *
    tauId *
    photonId *
    jetCorr *
    jetId *
    metCorr *
    preRun2Sync
)
