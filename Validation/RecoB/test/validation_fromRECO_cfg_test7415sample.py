#Basic example configration file to run the b-tagging validation sequence.
import FWCore.ParameterSet.Config as cms
process = cms.Process("validation")

"""
start customization
"""

#Enter here the Global tags
tag = '74X_mcRun2_asymptotic_v2'
#Do you want to apply JEC? For data, no need to add 'Residual', the code is checking if events are Data or MC and add 'Residual' for Data.
applyJEC = True
corrLabel = 'ak4PFCHSL1FastL2L3'
#Data or MC?
runOnMC    = True
#Flavour plots for MC: "all" = plots for all jets ; "dusg" = plots for d, u, s, dus, g independently ; not mandatory and any combinations are possible 
#b, c, light (dusg), non-identified (NI), PU jets plots are always produced
flavPlots = "allbcldusg"
#Check if jets originate from PU? option recommended (only for MC)
PUid = True
#List of taggers and taginfo to be considered (see example in: DQMOffline/RecoB/python/bTagCommon_cff.py)
from DQMOffline.RecoB.bTagCommon_cff import *
tagConfig = cms.VPSet(
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            folder = cms.string("CSVv2")
        ),
        #cms.PSet(
        #     bTagSimpleSVAnalysisBlock,
        #     label = cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"),
        #     folder = cms.string("SSVHE")
        #),

)

"""
end customization
"""

###prints###
print "is it MC ? : ", runOnMC
print "Global Tag : ", tag
############

process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.Core.DQM_cfg")
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

if runOnMC:
    #for MC jet flavour
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi")
    process.ak4JetFlavourInfos.jets = cms.InputTag("ak4PFJetsCHS")
    process.ak4JetFlavourInfos.hadronFlavourHasPriority = cms.bool(True)
    process.flavourSeq = cms.Sequence(
        process.selectedHadronsAndPartons *
        process.ak4JetFlavourInfos
    )
    #Validation sequence
    process.load("Validation.RecoB.bTagAnalysis_cfi")
    process.bTagValidation.jetMCSrc = 'ak4JetFlavourInfos'
    process.bTagValidation.tagConfig = tagConfig
    process.bTagHarvestMC.tagConfig = tagConfig
    process.bTagValidation.flavPlots = flavPlots
    process.bTagHarvestMC.flavPlots = flavPlots
    process.bTagValidation.doPUid = cms.bool(PUid)
    process.bTagValidation.doJEC = applyJEC
    process.bTagValidation.JECsource = cms.string(corrLabel)
    process.ak4GenJetsForPUid = cms.EDFilter("GenJetSelector",
                                             src = cms.InputTag("ak4GenJets"),
                                             cut = cms.string('pt > 8.'),
                                             filter = cms.bool(False)
                                             )
    process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
    process.patJetGenJetMatch.matched = cms.InputTag("ak4GenJetsForPUid")
    process.patJetGenJetMatch.maxDeltaR = cms.double(0.25)
    process.patJetGenJetMatch.resolveAmbiguities = cms.bool(True)
else :
    process.load("DQMOffline.RecoB.bTagAnalysisData_cfi")
    process.bTagAnalysis.tagConfig = tagConfig
    process.bTagHarvest.tagConfig = tagConfig
    process.bTagAnalysis.doJEC = applyJEC
    process.bTagAnalysis.JECsource = cms.string(corrLabel)

# load the full reconstraction configuration, to make sure we're getting all needed dependencies
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.GlobalTag.globaltag = tag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

if runOnMC:
    process.dqmSeq = cms.Sequence(process.ak4GenJetsForPUid * process.patJetGenJetMatch * process.flavourSeq * process.bTagValidation * process.bTagHarvestMC * process.dqmSaver)
else:
    process.dqmSeq = cms.Sequence(process.bTagAnalysis * process.bTagHarvest * process.dqmSaver)

process.plots = cms.Path(process.dqmSeq)
    
process.dqmEnv.subSystemFolder = 'BTAG'
process.dqmSaver.producer = 'DQM'
process.dqmSaver.workflow = '/POG/BTAG/BJET'
process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd =cms.untracked.bool(True) 
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.PoolSource.fileNames = [
       '/store/relval/CMSSW_7_4_12/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_asymptotic_v2-v1/00000/20A1016F-1D5B-E511-B4B9-00261894394D.root',
       '/store/relval/CMSSW_7_4_12/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_asymptotic_v2-v1/00000/28D2598A-205B-E511-915E-0025905A48D0.root',
       '/store/relval/CMSSW_7_4_12/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_asymptotic_v2-v1/00000/6E45DA91-205B-E511-A733-0025905A608C.root',
       '/store/relval/CMSSW_7_4_12/RelValTTbar_13/GEN-SIM-RECO/74X_mcRun2_asymptotic_v2-v1/00000/AA8A91CA-185B-E511-96AB-002618943810.root'
]

