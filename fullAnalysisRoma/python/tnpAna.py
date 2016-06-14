import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

isMC = False

process = cms.Process("tnpAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag import GlobalTag

if (isMC):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')     
elif (isMC==False):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v3', '')
print process.GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        # DY
        #"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring16_v1/diphotons_80_v1/DYToEE_NNPDF30_13TeV-powheg-pythia8/EXOSpring16_v1-diphotons_80_v1-v0-RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/160503_231305/0000/diphotonsMicroAOD_1.root"
        # data
        "/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring16_v1_p3/diphotons_80_v1/SingleElectron/EXOSpring16_v1_p3-diphotons_80_v1-v0-Run2016B-PromptReco-v2/160519_094902/0000/diphotonsMicroAOD_83.root"
        )
                            )

process.load("flashgg/MicroAOD/flashggPhotons_cfi")
process.load("flashgg/MicroAOD/flashggElectrons_cfi")
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TaP_output.root"))

process.tnpAna = cms.EDAnalyzer('TaPAnalyzer',
                                VertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                ElectronTag=cms.InputTag('flashggElectrons'),
                                genPhotonExtraTag = cms.InputTag("flashggGenPhotonsExtra"),    
                                DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                PileUpTag = cms.untracked.InputTag('slimmedAddPileupInfo'),
                                rhoTag = cms.InputTag('fixedGridRhoAll'),
                                rhoEleTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                bits = cms.InputTag("TriggerResults","","HLT"),
                                objects = cms.InputTag("selectedPatTrigger"),
                                MetTag=cms.InputTag('slimmedMETs'),
                                reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'), 
                                reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits'),
                                generatorInfo = cms.InputTag("generator"),
                                dopureweight = cms.untracked.int32(0),
                                sampleIndex  = cms.untracked.int32(10001), 
                                puWFileName  = cms.string('/afs/cern.ch/user/c/crovelli/public/json2015/prompt/singleEle/pileupWeights___processedAndGolden_2016B_june6__69mb.root'),
                                xsec         = cms.untracked.double(1),
                                kfac         = cms.untracked.double(1.),
                                sumDataset   = cms.untracked.double(1.)
                                )

process.p = cms.Path(process.tnpAna)

