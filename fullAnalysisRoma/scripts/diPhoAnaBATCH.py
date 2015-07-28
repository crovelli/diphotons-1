import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("diPhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = 'MCRUN2_74_V9A'         # 50ns 

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring( 

        # Spring15, RS 1TeV, k=0.1 
        "/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring15_v1/Spring15BetaV2/RSGravToGG_kMpl-01_M-1000_TuneCUEP8M1_13TeV-pythia8/EXOSpring15_v1-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150716_211524/0000/diphotonsMicroAOD_1.root" 
        )
                            )

process.load("flashgg/MicroAOD/flashggPhotons_cfi")
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")

process.TFileService = cms.Service("TFileService",fileName = cms.string("OUTPUT"))

process.diPhoAna = cms.EDAnalyzer('DiPhoAnalyzer',
                                  VertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                  genPhotonExtraTag = cms.InputTag("flashggGenPhotonsExtra"),   
                                  DiPhotonTag = cms.untracked.InputTag('flashggDiPhotons'),
                                  PileupTag = cms.untracked.InputTag('addPileupInfo'),
                                  generatorInfo = cms.InputTag("generator"),
                                  dopureweight = PU,
                                  sampleIndex  = SI,
                                  puWFileName  = weights,
                                  xsec         = XS,
                                  kfac         = KF,
                                  sumDataset   = SDS
                                  )

process.p = cms.Path(process.diPhoAna)

