import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag

from Configuration.Eras.Era_Phase2C11M9_cff import Phase2C11M9

process = cms.Process('MUNTUPLES',Phase2C11M9)

import subprocess
import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 '110X_mcRun4_realistic_v3',
                 #'110X_mcRun3_2021_realistic_v9',
                 #'auto:phase1_2021_realistic',
                 #'111X_upgrade2018_realistic_v1', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('nEvents',
                 #1000, #to run on a sub-sample
                 -1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Maximum number of processed events")

options.register('inputFile',                 'file:/store/relval/CMSSW_11_3_0_pre4/RelValZMM_14/GEN-SIM-RECO/113X_mcRun4_realistic_v4_2026D76noPU-v1/00000/08845437-06f6-4f0c-b267-912546e2cbbe.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "EOS folder with input files")

options.register('ntupleName',
                 'mutuple.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Folder and name ame for output ntuple")

options.parseArguments()

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

process.load('MuDPGAnalysis.MuonDPGNtuples.muNtupleProducer_cfi')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        numberOfThreads = cms.untracked.uint32(4))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))

#process.GlobalTag.globaltag = cms.string(options.globalTag)
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.source = cms.Source("PoolSource",                           
  fileNames = cms.untracked.vstring(options.inputFile),
  secondaryFileNames = cms.untracked.vstring()
)
 

process.TFileService = cms.Service('TFileService',
        fileName = cms.string(options.ntupleName)
    )



process.p = cms.Path(process.muNtupleProducer)


