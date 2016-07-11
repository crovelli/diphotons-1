import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("tnpAnaFlashgg")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_1_0-25ns_ICHEP16/2_1_0/SingleElectron/RunIISpring16DR80X-2_1_0-25ns_ICHEP16-2_1_0-v0-Run2016B-PromptReco-v1/160618_075833/0000/myMicroAODOutputFile_1.root" 
        )
                            )

from diphotons.MetaData.JobConfig import customize
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",1.)   ## chiara, era 1.e+3
customize.parse()

process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')     
if (customize.processType=="data"):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v3', '')
print process.GlobalTag


process.load("flashgg/MicroAOD/flashggPhotons_cfi")
process.load("flashgg/MicroAOD/flashggElectrons_cfi")
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TaP_output.root"))

process.tnpAna = cms.EDAnalyzer('TaPAnalyzer',
                                VertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                #ElectronTag=cms.InputTag('flashggElectrons'),
                                genPhotonExtraTag = cms.InputTag("flashggGenPhotonsExtra"),    
                                DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                PileUpTag = cms.untracked.InputTag('slimmedAddPileupInfo'),
                                rhoTag = cms.InputTag('fixedGridRhoAll'),
                                rhoEleTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                bits = cms.InputTag("TriggerResults","","HLT"),
                                objects = cms.InputTag("selectedPatTrigger"),
                                MetTag=cms.InputTag('slimmedMETs'),
                                generatorInfo = cms.InputTag("generator"),
                                dopureweight = cms.untracked.int32(1),
                                sampleIndex  = cms.untracked.int32(1),  
                                lumiWeight   = cms.untracked.double(1), 
                                puWFileName  = cms.string('/afs/cern.ch/user/c/crovelli/public/json2016/prompt/singleEle/pileupWeights___processedAndGolden_2016B_june22__69mb.root'),
                                )

process.p = cms.Path(process.tnpAna)

customize(process)
