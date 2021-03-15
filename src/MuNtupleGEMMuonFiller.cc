#include "MuDPGAnalysis/MuonDPGNtuples/src/MuNtupleGEMMuonFiller.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/Ref.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegment.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>

#include "TVectorF.h"
#include "TFile.h"

MuNtupleGEMMuonFiller::MuNtupleGEMMuonFiller(edm::ConsumesCollector && collector,
					     const std::shared_ptr<MuNtupleConfig> config, 
					     std::shared_ptr<TTree> tree, const std::string & label) : 
MuNtupleTrackBaseFiller(config, tree, label), m_nullVecF()
{

  edm::InputTag & muonTag = m_config->m_inputTags["muonTag"];
  if (muonTag.label() != "none") m_muToken = collector.consumes<reco::MuonCollection>(muonTag);

  edm::InputTag & primaryVerticesTag = m_config->m_inputTags["primaryVerticesTag"];
  if (primaryVerticesTag.label() != "none") m_primaryVerticesToken = collector.consumes<std::vector<reco::Vertex>>(primaryVerticesTag);

  edm::InputTag & gemSegmentTag = m_config->m_inputTags["gemSegmentTag"];
  if (gemSegmentTag.label() != "none") m_gemSegmentToken = collector.consumes<GEMSegmentCollection>(gemSegmentTag);
  
  edm::InputTag & cscSegmentTag = m_config->m_inputTags["cscSegmentTag"];
  if(cscSegmentTag.label() != "none") m_cscSegmentToken = collector.consumes<CSCSegmentCollection>(cscSegmentTag);

  edm::InputTag & trigResultsTag = m_config->m_inputTags["trigResultsTag"];
  if (trigResultsTag.label() != "none") m_trigResultsToken = collector.consumes<edm::TriggerResults>(trigResultsTag);

  edm::InputTag & trigEventTag = m_config->m_inputTags["trigEventTag"];
  if (trigEventTag.label() != "none") m_trigEventToken = collector.consumes<trigger::TriggerEvent>(trigEventTag);

  edm::InputTag & gemRecHitTag = m_config->m_inputTags["gemRecHitTag"];
  if (gemRecHitTag.label() != "none") m_gemRecHitToken = collector.consumes<GEMRecHitCollection>(gemRecHitTag);
  
}

MuNtupleGEMMuonFiller::~MuNtupleGEMMuonFiller() 
{ 

}

void MuNtupleGEMMuonFiller::initialize()
{

  m_tree->Branch((m_label + "_nMuons").c_str(), &m_nMuons);
  
  m_tree->Branch((m_label + "_pt").c_str(), &m_pt);
  m_tree->Branch((m_label + "_phi").c_str(), &m_phi);
  m_tree->Branch((m_label + "_eta").c_str(), &m_eta);
  m_tree->Branch((m_label + "_charge").c_str(), &m_charge);

  m_tree->Branch((m_label + "_isGlobal").c_str(), &m_isGlobal);
  m_tree->Branch((m_label + "_isStandalone").c_str(), &m_isStandalone);
  m_tree->Branch((m_label + "_isTracker").c_str(), &m_isTracker);
  m_tree->Branch((m_label + "_isGEM").c_str(), &m_isGEM);
  m_tree->Branch((m_label + "_isCSC").c_str(), &m_isCSC);
  m_tree->Branch((m_label + "_isME11").c_str(), &m_isME11);
  
  m_tree->Branch((m_label + "_isLoose").c_str(), &m_isLoose);
  m_tree->Branch((m_label + "_isMedium").c_str(), &m_isMedium);
  m_tree->Branch((m_label + "_isTight").c_str(), &m_isTight);

  m_tree->Branch((m_label + "_propagated_isME11").c_str(), &m_propagated_isME11);

  m_tree->Branch((m_label + "_path_length").c_str(), &m_path_length);

  m_tree->Branch((m_label + "_isinsideout").c_str(), &m_isinsideout);
  m_tree->Branch((m_label + "_isincoming").c_str(), &m_isincoming);

  m_tree->Branch((m_label + "_propagated_region").c_str(), &m_propagated_region);
  m_tree->Branch((m_label + "_propagated_layer").c_str(), &m_propagated_layer);
  m_tree->Branch((m_label + "_propagated_chamber").c_str(), &m_propagated_chamber);
  m_tree->Branch((m_label + "_propagated_etaP").c_str(), &m_propagated_etaP);

  m_tree->Branch((m_label + "_propagated_pt").c_str(), &m_propagated_pt);
  m_tree->Branch((m_label + "_propagated_phi").c_str(), &m_propagated_phi);
  m_tree->Branch((m_label + "_propagated_eta").c_str(), &m_propagated_eta);
  m_tree->Branch((m_label + "_propagated_charge").c_str(), &m_propagated_charge);

  m_tree->Branch((m_label + "_propagatedLoc_x").c_str(), &m_propagatedLoc_x);
  m_tree->Branch((m_label + "_propagatedLoc_y").c_str(), &m_propagatedLoc_y);
  m_tree->Branch((m_label + "_propagatedLoc_z").c_str(), &m_propagatedLoc_z);
  m_tree->Branch((m_label + "_propagatedLoc_r").c_str(), &m_propagatedLoc_r);
  m_tree->Branch((m_label + "_propagatedLoc_phi").c_str(), &m_propagatedLoc_phi);
  m_tree->Branch((m_label + "_propagatedLoc_errX").c_str(), &m_propagatedLoc_errX);
  m_tree->Branch((m_label + "_propagatedLoc_errY").c_str(), &m_propagatedLoc_errY);
  m_tree->Branch((m_label + "_propagatedGlb_x").c_str(), &m_propagatedGlb_x);
  m_tree->Branch((m_label + "_propagatedGlb_y").c_str(), &m_propagatedGlb_y);
  m_tree->Branch((m_label + "_propagatedGlb_z").c_str(), &m_propagatedGlb_z);
  m_tree->Branch((m_label + "_propagatedGlb_r").c_str(), &m_propagatedGlb_r);
  m_tree->Branch((m_label + "_propagatedGlb_phi").c_str(), &m_propagatedGlb_phi);
 
}

void MuNtupleGEMMuonFiller::clear()
{

  m_nMuons = 0;

  m_pt.clear();
  m_phi.clear();
  m_eta.clear();
  m_charge.clear();

  m_isGlobal.clear();
  m_isStandalone.clear();
  m_isTracker.clear();
  m_isGEM.clear();
  m_isCSC.clear();
  m_isME11.clear();
  
  m_isLoose.clear();
  m_isMedium.clear();
  m_isTight.clear();

  m_path_length = 0;
  
  m_isinsideout.clear();
  m_isincoming.clear();

  m_propagated_region.clear();
  m_propagated_layer.clear();
  m_propagated_chamber.clear();
  m_propagated_etaP.clear();
  m_propagated_isME11.clear();
  
  m_propagated_pt.clear();
  m_propagated_phi.clear();
  m_propagated_eta.clear();
  m_propagated_charge.clear();

  m_propagatedLoc_x.clear();
  m_propagatedLoc_y.clear();
  m_propagatedLoc_z.clear();
  m_propagatedLoc_r.clear();
  m_propagatedLoc_phi.clear();
  m_propagatedLoc_errX.clear();
  m_propagatedLoc_errY.clear();
  m_propagatedGlb_x.clear();
  m_propagatedGlb_y.clear();
  m_propagatedGlb_z.clear();
  m_propagatedGlb_r.clear();
  m_propagatedGlb_phi.clear();
      
}

void MuNtupleGEMMuonFiller::fill_new(const edm::Event & ev, const edm::EventSetup & environment)
{

  clear();

  auto muons = conditionalGet<reco::MuonCollection>(ev, m_muToken, "MuonCollection");
  auto gem_segments = conditionalGet<GEMSegmentCollection>(ev,m_gemSegmentToken, "GEMSegmentCollection");
  auto csc_segments = conditionalGet<CSCSegmentCollection>(ev,m_cscSegmentToken, "CSCSegmentCollection" );
  auto vtxs = conditionalGet<std::vector<reco::Vertex>>(ev, m_primaryVerticesToken, "std::vector<reco::Vertex>");

  bool isCSC = false;
  bool isME11 = false;
  
  edm::Handle<GEMRecHitCollection> rechit_collection;
  ev.getByToken(m_gemRecHitToken, rechit_collection);
  if (not rechit_collection.isValid()) {
    std::cout << "GEMRecHitCollection is invalid" << std::endl;
    return;
  }
 
  m_config->muon_service->update(environment);

  edm::ESHandle<Propagator>&& propagator_any = m_config->muon_service->propagator("SteppingHelixPropagatorAny");
  if (not propagator_any.isValid()) {
    std::cout<< "Any Propagator is invalid" << std::endl;
    return;
  }
  
  edm::ESHandle<Propagator>&& propagator_along = m_config->muon_service->propagator("SteppingHelixPropagatorAlong");
  if (not propagator_along.isValid()) {
    std::cout<< "Along Propagator is invalid" << std::endl;
    return;
  }

  edm::ESHandle<Propagator>&& propagator_opposite = m_config->muon_service->propagator("SteppingHelixPropagatorOpposite");
  if (not propagator_opposite.isValid()) {
    std::cout<< "Opposite Propagator is invalid" << std::endl;
    return;
  }

    
  edm::ESHandle<GEMGeometry> gem = m_config->m_gemGeometry;
  if (not gem.isValid()) {
    std::cout << "GEMGeometry is invalid" << std::endl;
    return;
  }

  edm::ESHandle<CSCGeometry> csc = m_config->m_cscGeometry;
  if (not csc.isValid()) {
    std::cout << "CSCGeometry is invalid" << std::endl;
    return;
  }
    
  edm::ESHandle<TransientTrackBuilder> transient_track_builder = m_config->m_transientTrackBuilder;
  if (not transient_track_builder.isValid()) {
    std::cout << "TransientTrack is invalid" << std::endl;
    return;
  }

  if (muons.isValid()) // && csc_segments.isValid() && vtxs.isValid()) 
    {
      //loop on recoMuons  
      for (const auto & muon : (*muons))
	{

	  m_pt.push_back(muon.pt());
	  m_eta.push_back(muon.eta());
	  m_phi.push_back(muon.phi());
	  m_charge.push_back(muon.charge());

	  m_isGlobal.push_back(muon.isGlobalMuon());
	  m_isStandalone.push_back(muon.isStandAloneMuon());
	  m_isTracker.push_back(muon.isTrackerMuon());
	  m_isGEM.push_back(muon.isGEMMuon());

	  m_isLoose.push_back(muon.passed(reco::Muon::CutBasedIdLoose));
	  m_isMedium.push_back(muon.passed(reco::Muon::CutBasedIdMedium));
	  m_isTight.push_back(muon.passed(reco::Muon::CutBasedIdTight));

	  m_nMuons++;

          isCSC = false;
          isME11 = false;

	  if(!muon.outerTrack().isNull())   //Cosmics
	  //if(!muon.globalTrack().isNull())   //Zmumu
	    {

	      //const reco::Track* track = muon.globalTrack().get();   //Zmumu
	      const reco::Track* track = muon.outerTrack().get();   //Cosmics
	      
	      if (track == nullptr) {
		std::cout << "failed to get muon track" << std::endl;
                continue;
              }
              
  
	      const reco::TrackRef outerTrackRef = muon.outerTrack();   //Cosmics
	      //const reco::TrackRef trackRef = muon.globalTrack();     //Zmumu
	    

	      float p2_in = track->innerMomentum().mag2();
	      float p2_out = track->outerMomentum().mag2();
	      float pos_out = track->outerPosition().mag2();
	      float pos_in = track->innerPosition().mag2();

	      bool is_insideout = pos_in > pos_out;

	      if(is_insideout)
		{
		  std::swap(pos_in, pos_out);
		  std::swap(p2_in, p2_out);
		}
	            
	      bool is_incoming = p2_out > p2_in;
	            
	      const reco::TransientTrack&& transient_track = transient_track_builder->build(track);
	      if (not transient_track.isValid()) 
		{
		  std::cout<<"failed  to build TransientTrack" << std::endl;
		  continue;
		}

	     	      
	      const auto&& start_state = is_insideout ? transient_track.outermostMeasurementState() : transient_track.innermostMeasurementState();
	      auto& propagator = is_incoming ? propagator_along : propagator_opposite;
	      
              //loop on recHits which form the outerTrack 
	      auto recHitMu = outerTrackRef->recHitsBegin();
	      auto recHitMuEnd = outerTrackRef->recHitsEnd();     //Cosmics

	      /*auto recHitMu = trackRef->recHitsBegin();
		auto recHitMuEnd = trackRef->recHitsEnd();*/
	      
	      for(; recHitMu != recHitMuEnd; ++recHitMu)
		{
		  DetId detId = (*recHitMu)->geographicalId();
		    
		  if(detId.det() == DetId::Muon && detId.subdetId() == MuonSubdetId::CSC)
		    {
		      isCSC = true;
		      
		      const CSCDetId csc_id{detId};
		      
		      /*std::cout << "rec csc evento " << std::endl;
			      const CSCRecHit2D *cscRecHitSta = dynamic_cast<const CSCRecHit2D *>(*recHitMu);
			      if (cscRecHitSta == nullptr)
				    {
				    std::cout << "dynamic cast failure" << std::endl;
				      continue;
				      }
		      */

		      /*const CSCSegment *cscSegmentSta = dynamic_cast<const CSCSegment*>(*recHitMu);
			      if(cscSegmentSta == nullptr)
			      {
			      std::cout << "dynamic cast segment failure" << std::endl;
			        //continue;
				}
		      */
		      if(csc_id.station() == 1 && ((csc_id.ring() == 1) || (csc_id.ring() == 4)) ) isME11 = true;
		    }
		} //loop on recHits
	      m_isCSC.push_back(isCSC);
	      m_isME11.push_back(isME11);
	           
              //if at least one CSC hit is found, perform propagation 
	      if(isCSC)
                { 
                  for (const GEMRegion* gem_region : gem->regions())
                    {
                      bool is_opposite_region = muon.eta() * gem_region->region() < 0;
                      if (is_incoming xor is_opposite_region) continue;
                      
                      for (const GEMStation* station : gem_region->stations())
                        {
			  //std::cout << "station" << std::endl;
                          for (const GEMSuperChamber* super_chamber : station->superChambers())
                            {
                              for (const GEMChamber* chamber : super_chamber->chambers())
                                {
				  
                                  const BoundPlane& bound_plane = chamber->surface();
                                        
                                  const auto& dest_state = propagator_any->propagate(start_state, bound_plane);
                                  if (not dest_state.isValid())
                                    {
				      std::cout << "failed to propagate" << std::endl;
                                      continue;
                                    }

                                  const GlobalPoint&& dest_global_pos = dest_state.globalPosition();
                                  
                                  const GEMEtaPartition* eta_partition = findEtaPartition(chamber, dest_global_pos);
                                  if (eta_partition == nullptr)
                                    {
				      //std::cout << "failed to find GEMEtaPartition" << std::endl; 
                                      continue;
                                    }

                                  const GEMDetId&& gem_id = eta_partition->id();

                                  m_propagated_isME11.push_back(isME11);

                                  const LocalPoint&& dest_local_pos = eta_partition->toLocal(dest_global_pos);
				  const LocalError&& dest_local_err = dest_state.localError().positionError();

                                  m_propagatedGlb_x.push_back(dest_global_pos.x());
                                  m_propagatedGlb_y.push_back(dest_global_pos.y());
                                  m_propagatedGlb_z.push_back(dest_global_pos.z());
                                  m_propagatedGlb_r.push_back(dest_global_pos.perp());
                                  m_propagatedGlb_phi.push_back(dest_global_pos.phi());
                                  
                                  m_propagated_pt.push_back(muon.pt());
                                  m_propagated_phi.push_back(muon.phi());
                                  m_propagated_eta.push_back(muon.eta());
                                  m_propagated_charge.push_back(muon.charge());
                        
                                  m_propagatedLoc_x.push_back(dest_local_pos.x());
                                  m_propagatedLoc_phi.push_back(dest_local_pos.phi());
                                  m_propagatedLoc_r.push_back(dest_local_pos.perp());
                                  
				  m_propagatedLoc_errX.push_back(dest_local_err.xx());
				  m_propagatedLoc_errY.push_back(dest_local_err.yy());
				  
                                  m_propagated_region.push_back(gem_id.region());
                                  m_propagated_layer.push_back(gem_id.layer());
                                  m_propagated_chamber.push_back(gem_id.chamber());
                                  m_propagated_etaP.push_back(gem_id.roll());
                                  
                                  m_isinsideout.push_back(is_insideout);
                                  m_isincoming.push_back(is_incoming);
				}  
			    }
			}   
		    }
		}//isCSC
	    }//!muon.outerTrack().isNull()
	}//loop on reco muons
    }//muons.isValid()
  
  return;
  
}


const GEMEtaPartition* MuNtupleGEMMuonFiller::findEtaPartition(const GEMChamber* chamber, const GlobalPoint& global_point) {

  for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
    const LocalPoint&& local_point = eta_partition->toLocal(global_point);
    const LocalPoint local_point_2d(local_point.x(), local_point.y(), 0.0f);
    if (eta_partition->surface().bounds().inside(local_point_2d)) 
      return eta_partition;
  }

  return nullptr;
}
