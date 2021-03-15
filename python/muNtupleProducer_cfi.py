import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy

muNtupleProducer = cms.EDAnalyzer("MuNtupleProducer",
                                  MuonServiceProxy,
                                  residualXCut = cms.double(5.0),
                                  ph1DtDigiTag = cms.untracked.InputTag("muonDTDigis"),
                                  ph2DtDigiTag = cms.untracked.InputTag("none"),

                                  ph1DtSegmentTag = cms.untracked.InputTag("dt4DSegments"),        
                                  ph2DtSegmentTag = cms.untracked.InputTag("none"),

                                  ph1DTtTrigMode = cms.untracked.string('DTTTrigSyncFromDB'),
                                  gemDigiTag = cms.untracked.InputTag("muonGEMDigis"),
                                  gemRecHitTag = cms.untracked.InputTag("gemRecHits"),
                                  gemSegmentTag = cms.untracked.InputTag("gemSegments"),
                                  cscSegmentTag = cms.untracked.InputTag("cscSegments"),
                                  muonTag = cms.untracked.InputTag("muons"),
                                  primaryVerticesTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                                  ph1DTtTrigModeConfig = cms.untracked.PSet(vPropWire = cms.double(24.4),
                                                                            doTOFCorrection = cms.bool(False),
                                                                            tofCorrType = cms.int32(2),
                                                                            wirePropCorrType = cms.int32(0),
                                                                            doWirePropCorrection = cms.bool(False),
                                                                            doT0Correction = cms.bool(True),
                                                                            tTrigLabel = cms.string(''),
                                                                            t0Label = cms.string(''),
                                                                            debug = cms.untracked.bool(False)
                                                                        )
)

