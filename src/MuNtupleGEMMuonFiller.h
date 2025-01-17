#ifndef MuNtuple_MuNtupleGEMMuonFiller_h
#define MuNtuple_MuNtupleGEMMuonFiller_h

#include "MuDPGAnalysis/MuonDPGNtuples/src/MuNtupleBaseFiller.h"
#include "MuDPGAnalysis/MuonDPGNtuples/src/MuNtupleTrackBaseFiller.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/GEMRecHit/interface/GEMSegment.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include <vector>
#include <fstream>
#include "TClonesArray.h"
#include "TVectorF.h"
#include "TFile.h"

class MuNtupleGEMMuonFiller : public MuNtupleTrackBaseFiller
{

 public:

  /// Constructor
  MuNtupleGEMMuonFiller(edm::ConsumesCollector && collector,
		     const std::shared_ptr<MuNtupleConfig> config, 
		     std::shared_ptr<TTree> tree, const std::string & label);
  
  ///Destructor
  virtual ~MuNtupleGEMMuonFiller();
  
  /// Intialize function : setup tree branches etc ... 
  virtual void initialize() final;
  
  /// Clear branches before event filling 
  virtual void clear() final;
 
  virtual void fill_new(const edm::Event & ev, const edm::EventSetup & environment) final;
  
 private:

  edm::EDGetTokenT<reco::MuonCollection>      m_muToken;
  edm::EDGetTokenT<GEMSegmentCollection>  m_gemSegmentToken;
  edm::EDGetTokenT<std::vector<reco::Vertex>> m_primaryVerticesToken;
  edm::EDGetTokenT<CSCSegmentCollection> m_cscSegmentToken;
  edm::EDGetTokenT<GEMRecHitCollection> m_gemRecHitToken;

  edm::EDGetTokenT<edm::TriggerResults>   m_trigResultsToken;
  edm::EDGetTokenT<trigger::TriggerEvent> m_trigEventToken;

  const GEMRecHit *findMatchedHit(const float,  const GEMRecHitCollection::range );
  const GEMEtaPartition*  findEtaPartition(const GEMChamber*, const GlobalPoint&);

  const TrackingRecHit *getHitPtr(edm::OwnVector<TrackingRecHit>::const_iterator iter) const {return &*iter; }
  const TrackingRecHit *getHitPtr(const trackingRecHit_iterator &iter) const {return &**iter; }

  TVectorF m_nullVecF;

  std::fstream outdata;
  
  unsigned int m_nMuons; 

  std::vector<float> m_pt;     // muon pT [GeV/c]
  std::vector<float> m_phi;    // muon phi [rad]
  std::vector<float> m_eta;    // muon eta
  std::vector<short> m_charge; // muon charge

  std::vector<bool>  m_isGlobal;
  std::vector<bool>  m_isStandalone;
  std::vector<bool>  m_isTracker;
  //std::vector<bool>  m_isTrackerArb;
  std::vector<bool>  m_isGEM;
  std::vector<bool>  m_isCSC;
  std::vector<bool> m_isME11;
  
  std::vector<bool> m_propagated_isME11;

  std::vector<bool>  m_isLoose;  // Loose muon ID
  std::vector<bool>  m_isMedium; // Medium muon ID
  std::vector<bool>  m_isTight;  // Tight muon ID

  std::vector<bool> m_isincoming;
  std::vector<bool> m_isinsideout;
  
  float m_path_length;

  std::vector<int> m_propagated_region;
  std::vector<int> m_propagated_layer;
  std::vector<int> m_propagated_chamber;
  std::vector<int> m_propagated_etaP;

  std::vector<float> m_propagated_pt;
  std::vector<float> m_propagated_phi;
  std::vector<float> m_propagated_eta;
  std::vector<float> m_propagated_charge;

  std::vector<float> m_propagatedLoc_x;
  std::vector<float> m_propagatedLoc_errX;
  std::vector<float> m_propagatedLoc_y;
  std::vector<float> m_propagatedLoc_errY;
  std::vector<float> m_propagatedLoc_z;
  std::vector<float> m_propagatedLoc_r;
  std::vector<float> m_propagatedLoc_phi;
  std::vector<float> m_propagatedGlb_x;
  std::vector<float> m_propagatedGlb_y;
  std::vector<float> m_propagatedGlb_z;
  std::vector<float> m_propagatedGlb_r;
  std::vector<float> m_propagatedGlb_phi;
  
};


#endif
