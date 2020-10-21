import FWCore.ParameterSet.Config as cms
import copy

# Define collections
muon             = 'slimmedMuons'
electron        = 'slimmedElectrons'

MC=cms.bool(True)
MC_Signal=cms.bool(True)
DATA=cms.bool(False)

from Configuration.StandardSequences.Eras import eras

#process = cms.Process("Demo")
#process = cms.Process("CTPPSTestProtonReconstruction", eras.ctpps_2016)
# load common code
from direct_simu_reco_cff import *
process = cms.Process('CTPPSTestAcceptance', era)
process.load("direct_simu_reco_cff")
SetDefaults(process)



process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(35000) )

process.options=cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:signal_miniAOD.root"),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
if MC:
	process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3', '') # MC
else:
	process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v11') # data


#################################
  ###  JET TOOLBOX FOR CHS ###
#################################
'''
from Configuration.Eras.Modifier_run2_muon_2016_cff import run2_muon_2016
from Configuration.Eras.Modifier_run2_egamma_2017_cff import run2_egamma_2017
from Configuration.Eras.Modifier_run2_egamma_2016_cff import run2_egamma_2016
from Configuration.Eras.Modifier_run2_tau_ul_2016_cff import run2_tau_ul_2016
from Configuration.Eras.Modifier_ctpps_2017_cff import ctpps_2017
from Configuration.Eras.Modifier_pixel_2016_cff import pixel_2016
'''
# AK R=0.8 jets from CHS inputs with basic grooming, W tagging, and top tagging                                                            
from JMEAnalysis.JetToolbox.jetToolbox_cff import *
#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox( process, 'ak8', 'ak8JetSubs', 'noOutput',
#jetToolbox( process, 'ak8', 'jetSequence', 'noOutput',
                PUMethod='CHS',runOnMC=MC,
                addPruning=True, addSoftDrop=False ,           # add basic grooming                                                            
                addTrimming=False, addFiltering=False,
                addSoftDropSubjets=False,
                addNsub=True, maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4                                      
                Cut='pt > 100.0',
                bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
                #bTagDiscriminators=['pfCombinedSecondaryVertexV2BJetTags'],
                # added L1FastJet on top of the example config file
                JETCorrPayload = 'AK8PFchs', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
                )
jetToolbox( process, 'ak4', 'ak4JetSubs', 'noOutput',
#jetToolbox( process, 'ak4', 'jetSequence', 'noOutput',
                PUMethod='CHS',runOnMC=MC,
                addPruning=True, addSoftDrop=False ,           # add basic grooming                                                            
                addTrimming=False, addFiltering=False,
                addSoftDropSubjets=False,
                addNsub=True, maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4                                      
                Cut='pt > 10.0',
                bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
                #bTagDiscriminators=['pfCombinedSecondaryVertexV2BJetTags'],
                # added L1FastJet on top of the example config file
                JETCorrPayload = 'AK4PFchs', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
                )





if MC:
	JetAK4Collection="selectedPatJetsAK4PFCHS" #"slimmedJets"
	JetAK8Collection="selectedPatJetsAK8PFCHS" #"slimmedJetsAK8"

###################    C L E A N E R    ###################
process.cleanPatMuons = cms.EDProducer("PATMuonCleaner",
    src = cms.InputTag(muon),
    # preselection (any string-based cut for pat::Muon)
    preselection = cms.string("passed('CutBasedIdTight') && passed('PFIsoTight') && pt() > 50 && abs(eta()) < 2.4"),
#    preselection = cms.string("muonID('GlobalMuonPromptTight') && numberOfMatchedStations() > 1 && dB() < 0.2 && pt() > 50 && abs(eta()) < 2.4 && ((pfIsolationR04().sumChargedHadronPt + max(0., pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt))/pt() < 0.1)"),  # "&& ((isolationR03().sumPt)/(pt()) < 0.1)",
    # overlap checking configurables
    checkOverlaps = cms.PSet(),
    # finalCut (any string-based cut for pat::Muon)
    finalCut = cms.string(''),
)

process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
    ## pat electron input source
    src = cms.InputTag(electron),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string("electronID('cutBasedElectronID-Summer16-80X-V1-tight') && pt() > 50 && abs(eta()) < 2.4"), #electronID('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight
        #(pfIsolationVariables().sumChargedHadronPt()+max(0.0,pfIsolationVariables().sumNeutralHadronEt()+pfIsolationVariables().sumPhotonEt()-rho*effArea))/pt())
    # overlap checking configurables
    checkOverlaps = cms.PSet(),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)

process.cleanJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag(JetAK4Collection),

    # preselection (any string-based cut on pat::Jet)
    preselection = cms.string(''),

    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("cleanPatMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
        electrons = cms.PSet(
           src       = cms.InputTag("cleanPatElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
    ),
    # finalCut (any string-based cut on pat::Jet)
    finalCut = cms.string(''),
)
process.cleanJetsAK8 = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag(JetAK8Collection), # selectedPatJetsAK8PFCHS
    # preselection (any string-based cut on pat::Jet)
    preselection = cms.string(''),

    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("cleanPatMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(1.0),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
        electrons = cms.PSet(
           src       = cms.InputTag("cleanPatElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(1.0),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
    ),
    # finalCut (any string-based cut on pat::Jet)
    finalCut = cms.string(''),
)
###########################################################


###################  J E T  F I L T E R   #################
process.filterJets = cms.EDFilter( 'FilterAK8Jet'
        , jetsAK8                  = cms.InputTag(JetAK8Collection)
)
###########################################################


# update settings of beam-smearing module
process.beamDivergenceVtxGenerator.src = cms.InputTag("")
process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
#  cms.InputTag("genPUProtons", ""),
  cms.InputTag("prunedGenParticles")
#  cms.InputTag("genParticles")
)

# do not apply vertex smearing again
process.ctppsBeamParametersESSource.vtxStddevX = 0
process.ctppsBeamParametersESSource.vtxStddevY = 0
process.ctppsBeamParametersESSource.vtxStddevZ = 0

# undo CMS vertex shift
process.ctppsBeamParametersESSource.vtxOffsetX45 = -1.048 * 1E-1
process.ctppsBeamParametersESSource.vtxOffsetY45 = -1.687 * 1E-1
process.ctppsBeamParametersESSource.vtxOffsetZ45 = +12.04 * 1E-1
#process.ctppsBeamParametersESSource.vtxOffsetX45 = +0.2475 * 1E-1
#process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.6924 * 1E-1
#process.ctppsBeamParametersESSource.vtxOffsetZ45 = -8.1100 * 1E-1

process.demo = cms.EDAnalyzer('MakeNTuple'
        , jetsAK8                  = cms.InputTag('cleanJetsAK8')
        , jetsAK4                  = cms.InputTag('cleanJets')
	, genJets		   = cms.InputTag('slimmedGenJets')
        , genJetsAK8               = cms.InputTag('slimmedGenJetsAK8')
        , muons                    = cms.InputTag(muon)
        , electrons                = cms.InputTag(electron)
        , MET                      = cms.InputTag("slimmedMETs")
        , PFCand                   = cms.InputTag('packedPFCandidates')
        , PileupSumInfoInputTag = cms.InputTag('slimmedAddPileupInfo')
        , vertices                      = cms.InputTag('offlineSlimmedPrimaryVertices')
	, MC			   = MC
        , MC_Signal                = MC_Signal
        , DATA                     = DATA
        , ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP")
        , ppsRecoProtonMultiRPTag  = cms.InputTag("ctppsProtons", "multiRP")
	, TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')
)


process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string("out.root")
                                                                        )

process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(False)

if MC:
#	process.p = cms.Path(process.ecalBadCalibReducedMINIAODFilter*process.slimmedAK8JetsSmeared*process.slimmedAK4JetsSmeared*process.cleanPatElectrons*process.cleanPatMuons*process.cleanJets*process.cleanJetsAK8*process.genMet*process.uncorrectedMet*process.uncorrectedPatMet*process.Type1CorrForNewJEC*process.slimmedMETsNewJEC*process.shiftedMETCorrModuleForSmearedJets*process.slimmedMETsSmeared*process.filterJets*process.beamDivergenceVtxGenerator*process.ctppsDirectProtonSimulation*process.reco_local*process.ctppsProtons*process.demo)
        process.p = cms.Path(process.cleanPatElectrons*process.cleanPatMuons*process.cleanJets*process.cleanJetsAK8*process.filterJets*process.beamDivergenceVtxGenerator*process.ctppsDirectProtonSimulation*process.reco_local*process.ctppsProtons*process.demo)
else:
        process.p = cms.Path(process.hltFilter*process.METFilter*process.ecalBadCalibReducedMINIAODFilter*process.cleanPatElectrons*process.cleanPatMuons*process.cleanJets*process.cleanJetsAK8*process.uncorrectedMet*process.Type1CorrForNewJEC*process.slimmedMETsNewJEC*process.filterJets*process.totemRPUVPatternFinder*process.totemRPLocalTrackFitter*process.ctppsLocalTrackLiteProducer*process.ctppsProtons*process.demo)






