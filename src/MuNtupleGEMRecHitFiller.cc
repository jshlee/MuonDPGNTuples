#include "MuDPGAnalysis/MuonDPGNtuples/src/MuNtupleGEMRecHitFiller.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"


MuNtupleGEMRecHitFiller::MuNtupleGEMRecHitFiller(edm::ConsumesCollector && collector,
                                                 const std::shared_ptr<MuNtupleConfig> config,
                                                 std::shared_ptr<TTree> tree, const std::string & label) :

  MuNtupleBaseFiller(config, tree, label)
{

  edm::InputTag &iTag = m_config->m_inputTags["gemRecHitTag"];
  if (iTag.label() != "none") m_rechit_token_ = collector.consumes<GEMRecHitCollection>(iTag);

}

MuNtupleGEMRecHitFiller::~MuNtupleGEMRecHitFiller()
{

};

void MuNtupleGEMRecHitFiller::initialize()
{
  
  m_tree->Branch((m_label + "_nRecHits").c_str(), &m_nRecHits, (m_label + "_nRecHits/i").c_str());

  m_tree->Branch((m_label + "_cluster_size").c_str(), &m_rechit_cluster_size);
  m_tree->Branch((m_label + "_bx").c_str(), &m_rechit_bx);

  m_tree->Branch((m_label + "_region").c_str(), &m_rechit_region);
  m_tree->Branch((m_label + "_station").c_str(), &m_rechit_station);
  m_tree->Branch((m_label + "_chamber").c_str(), &m_rechit_chamber);
  m_tree->Branch((m_label + "_layer").c_str(), &m_rechit_layer);
  m_tree->Branch((m_label + "_etaPartition").c_str(), &m_rechit_etaPartition);

  m_tree->Branch((m_label + "_loc_r").c_str(), &m_rechit_loc_r);
  m_tree->Branch((m_label + "_loc_phi").c_str(), &m_rechit_loc_phi);
  m_tree->Branch((m_label + "_loc_x").c_str(), &m_rechit_loc_x);
  m_tree->Branch((m_label + "_loc_y").c_str(), &m_rechit_loc_y);
  m_tree->Branch((m_label + "_loc_z").c_str(), &m_rechit_loc_z);

  m_tree->Branch((m_label + "_g_r").c_str(), &m_rechit_g_r);
  m_tree->Branch((m_label + "_g_phi").c_str(), &m_rechit_g_phi);
  m_tree->Branch((m_label + "_g_x").c_str(), &m_rechit_g_x);
  m_tree->Branch((m_label + "_g_y").c_str(), &m_rechit_g_y);
  m_tree->Branch((m_label + "_g_z").c_str(), &m_rechit_g_z);
 
}

void MuNtupleGEMRecHitFiller::clear()
{

  m_nRecHits = 0;

  m_rechit_cluster_size.clear();
  m_rechit_bx.clear();

  m_rechit_region.clear();
  m_rechit_station.clear();
  m_rechit_chamber.clear();
  m_rechit_layer.clear();
  m_rechit_etaPartition.clear();

  m_rechit_loc_r.clear();
  m_rechit_loc_phi.clear();
  m_rechit_loc_x.clear();
  m_rechit_loc_y.clear();
  m_rechit_loc_z.clear();
 
  m_rechit_g_r.clear();
  m_rechit_g_phi.clear();
  m_rechit_g_x.clear();
  m_rechit_g_y.clear();
  m_rechit_g_z.clear();
  
  
}

void MuNtupleGEMRecHitFiller::fill(const edm::Event & ev)
{
  
  clear();

  edm::Handle<GEMRecHitCollection> rechit_collection;
  ev.getByToken(m_rechit_token_,rechit_collection);
  const auto gem = m_config->m_gemGeometry; 

    
  if (rechit_collection.isValid())
    {
	
      for (auto rechit = rechit_collection->begin(); rechit != rechit_collection->end(); rechit++)
        {
          GEMDetId gem_id{rechit->gemId()};
          // std::cout << "gemRechit" << gem_id << std::endl;
          m_rechit_etaPartition.push_back(gem_id.roll());
          m_rechit_region.push_back(gem_id.region());
          m_rechit_station.push_back(gem_id.station());
          m_rechit_chamber.push_back(gem_id.chamber());
          m_rechit_layer.push_back(gem_id.layer());

          const BoundPlane& surface = gem->idToDet(gem_id)->surface();
          GlobalPoint&& rechit_global_pos = surface.toGlobal(rechit->localPosition());
          //LocalPoint&& rechit_local_pos = surface.toLocal(rechit->globalPosition());
	    
          m_rechit_loc_r.push_back(rechit->localPosition().perp());
          m_rechit_loc_phi.push_back(rechit->localPosition().phi());
          m_rechit_loc_x.push_back(rechit->localPosition().x());
          m_rechit_loc_y.push_back(rechit->localPosition().y());
          m_rechit_loc_z.push_back(rechit->localPosition().z());

          m_rechit_g_r.push_back(rechit_global_pos.perp());
          m_rechit_g_phi.push_back(rechit_global_pos.phi());
          m_rechit_g_x.push_back(rechit_global_pos.x());
          m_rechit_g_y.push_back(rechit_global_pos.y());
          m_rechit_g_z.push_back(rechit_global_pos.z());

          //std::cout << "rec_x" << rechit->localPosition().x() << std::endl;
          //std::cout << "rec_y" << rechit->localPosition().y() << std::endl;
	    
          auto cluster_size = rechit->clusterSize();
          int bx = rechit->BunchX();
	    
          m_rechit_cluster_size.push_back(cluster_size);
          m_rechit_bx.push_back(bx);
          //std::cout << " cluster size: " << cluster_size << std::endl;

          m_nRecHits++;
        }
      //std::cout << " # rechits: " << m_nRecHits << std::endl;

    }

  return;

}
