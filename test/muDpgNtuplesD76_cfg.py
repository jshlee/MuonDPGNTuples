import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11M9_cff import Phase2C11M9
process = cms.Process('MUNTUPLES',Phase2C11M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration/StandardSequences/GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')

process.GlobalTag.globaltag = '113X_mcRun4_realistic_v4'
print 'process.GlobalTag=', process.GlobalTag

#process.options   = cms.untracked.PSet()
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",                           
    fileNames = cms.untracked.vstring('file:/store/relval/CMSSW_11_3_0_pre4/RelValZMM_14/GEN-SIM-RECO/113X_mcRun4_realistic_v4_2026D76noPU-v1/00000/08845437-06f6-4f0c-b267-912546e2cbbe.root')
)
 
process.TFileService = cms.Service('TFileService', fileName = cms.string('mutuple.root') )

process.load('MuDPGAnalysis.MuonDPGNtuples.muNtupleProducer_cfi')
process.p = cms.Path(process.muNtupleProducer)

