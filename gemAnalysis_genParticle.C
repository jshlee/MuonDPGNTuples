#define gemAnalysis_genParticle_cxx
#include "gemAnalysis_genParticle.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TCut.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <string.h>


gemAnalysis_genParticle::gemAnalysis_genParticle(const TString & inFileName,
                         const TString & outFileName) :
  m_inFile(inFileName,"READ"),m_outFile(outFileName,"RECREATE"),fChain(0)
{

  fChain = static_cast<TTree*>(m_inFile.Get("muNtupleProducer/MuDPGTree"));
  Init(fChain);

}


void gemAnalysis_genParticle::Loop()
{

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  book();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      
    std::cout << "event" << std::endl;
   
    for(std::size_t iMu = 0; iMu < mu_propagatedLoc_x->size(); ++iMu) //local matching
	    {

	      bool same_pt, same_eta, same_phi, matching_genprop;
	      float min_difference_pt_rel = 0.1;
	      float min_difference_dR = 0.05;
	      float diff_pt_rel = 0., diff_eta = 0., diff_phi = 0., diff_R = 0.;

	      //bool eta_cut = TMath::Abs(mu_propagated_eta->at(iMu)) >= 1.55 && TMath::Abs(mu_propagated_eta->at(iMu)) <= 2.18;
	      bool pt_cut = mu_propagated_pt->at(iMu) >= 2;
	      std::string iEta = "";
	      
	      same_pt = false;
	      same_eta = false;
	      same_phi = false;
	      matching_genprop = false;
	      
	      for(std::size_t iGen = 0; iGen < genParticle_PdgId->size(); ++iGen )
		{
		  bool pt_cut_gen = genParticle_Pt->at(iGen) >= 2;

		  if(pt_cut_gen && pt_cut && TMath::Abs(genParticle_PdgId->at(iGen)) == 13 && genParticle_MotherPdgId->at(iGen) == 23)
		    {
		      diff_pt_rel = (mu_propagated_pt->at(iMu) - genParticle_Pt->at(iGen))/mu_propagated_pt->at(iMu);
		      //diff_eta = mu_propagated_eta->at(iMu) - genParticle_Eta->at(iGen);
		      //diff_phi = mu_propagated_phi->at(iMu) - genParticle_Phi->at(iGen);
		      
		      diff_R = dR(mu_propagated_eta->at(iMu),genParticle_Eta->at(iGen),mu_propagated_phi->at(iMu),genParticle_Phi->at(iGen));

		      if(TMath::Abs(diff_pt_rel) < min_difference_pt_rel && diff_R < min_difference_dR)
			{
			  min_difference_pt_rel = diff_pt_rel;
			  min_difference_dR = diff_R;
			  matching_genprop = true;
			}
		      else continue;
		    }
		}
		  
	      if(matching_genprop && mu_propagated_region->at(iMu) == -1  && mu_propagated_isME11->at(iMu) && mu_propagated_pt->at(iMu) >= 2)
		    {
		      iEta = std::to_string(mu_propagated_etaP->at(iMu));
		      m_plots["MuPtPropagated_endcapneg_iEta_"+iEta]->Fill(mu_propagated_pt->at(iMu));
		      //m_plots["MuEtaPropagated_endcapneg_iEta_"+iEta]->Fill(mu_propagated_eta->at(iMu));
		      //m_plots["MuPhiPropagated_endcapneg_iEta_"+iEta]->Fill(mu_propagated_phi->at(iMu));


		      int j=0; //scanning index

                      for(j=0; j<10; j++)
                        {

			  float min_residual_x = 2 - (float)j*2/10;
			  std::string string = std::to_string(min_residual_x);
			  bool matching = false;
			  
			  float min_recHitPosition = 0;
			  float mupt = 0;
			  float mueta = 0;
			  float muphi = 0;
			  int etaP = 0;
			  bool isME11 = 0;
			  std::string iEtaRechit = "";
			  
			  for(std::size_t iRecHit = 0; iRecHit < gemRecHit_nRecHits; iRecHit++)
			    {
			      float delta_x = mu_propagatedLoc_x->at(iMu)-gemRecHit_loc_x->at(iRecHit);
			      bool same_eta = (mu_propagated_etaP->at(iMu) == gemRecHit_etaPartition->at(iRecHit));
			      bool same_chamber = (mu_propagated_chamber->at(iMu) ==  gemRecHit_chamber->at(iRecHit));
			      bool same_layer = ( mu_propagated_layer->at(iMu) == gemRecHit_layer->at(iRecHit) );
			      
			      if(std::fabs(delta_x) < min_residual_x &&  gemRecHit_region->at(iRecHit) == -1 && same_eta && same_chamber && same_layer)
				{
				  min_residual_x = delta_x;
				  min_recHitPosition = gemRecHit_loc_x->at(iRecHit);
				  etaP = gemRecHit_etaPartition->at(iRecHit);
				  iEtaRechit = std::to_string(etaP);
				  mupt = mu_propagated_pt->at(iMu);
				  mueta = TMath::Abs(mu_propagated_eta->at(iMu));
				  muphi = mu_propagated_phi->at(iMu);
				  
				  matching = true;
				}
			      else continue;
			    }//end loop rechit
			  
			  if(matching == true)
			    {
			      m_plots["ResidualPlotLocalX_endcapneg_dX_"+string]->Fill(min_residual_x);
			      m_plots["MuPtMatchedLocal_endcapneg_dX_"+string+"_iEta_"+iEta]->Fill(mupt);
			      //m_plots["MuEtaMatchedLocal_endcapneg_dX_"+string+"_iEta_"+iEta]->Fill(mueta);
			      //m_plots["MuPhiMatchedLocal_endcapneg_dX_"+string+"_iEta_"+iEta]->Fill(muphi);
			      
			    }
			  else continue;
			}//end scanning iEta
		    
		    }// if negative endcap 
	      
	      else if(matching_genprop && mu_propagated_region->at(iMu) == 1  && mu_propagated_isME11->at(iMu) && mu_propagated_pt->at(iMu) >= 2 )
		{
		  iEta = std::to_string(mu_propagated_etaP->at(iMu));
		  //m_plots["MuPtPropagated_endcappos_iEta_"+iEta]->Fill(mu_propagated_pt->at(iMu));
		  //m_plots["MuEtaPropagated_endcappos_iEta_"+iEta]->Fill(mu_propagated_eta->at(iMu));
		  //m_plots["MuPhiPropagated_endcapos_iEta_"+iEta]->Fill(mu_propagated_phi->at(iMu));


		  int j=0; //scanning index

		  for(j=0; j<10; j++)
		    {

		      float min_residual_x = 2 - (float)j*2/10;
		      std::string string = std::to_string(min_residual_x);
		      bool matching = false;

		      float min_recHitPosition = 0;
		      float mupt = 0;
		      float mueta = 0;
		      float muphi = 0;
		      int etaP = 0;
		      bool isME11 = 0;
		      std::string iEtaRechit = "";

		      for(std::size_t iRecHit = 0; iRecHit < gemRecHit_nRecHits; iRecHit++)
			{
			  float delta_x = mu_propagatedLoc_x->at(iMu)-gemRecHit_loc_x->at(iRecHit);
			  bool same_eta = (mu_propagated_etaP->at(iMu) == gemRecHit_etaPartition->at(iRecHit));
			  bool same_chamber = (mu_propagated_chamber->at(iMu) ==  gemRecHit_chamber->at(iRecHit));
			  bool same_layer = ( mu_propagated_layer->at(iMu) == gemRecHit_layer->at(iRecHit) );

			  if(std::fabs(delta_x) < min_residual_x &&  gemRecHit_region->at(iRecHit) == -1 && same_eta && same_chamber && same_layer)// && mu_propagated_isME11->at(iMu))
			    {
			      min_residual_x = delta_x;
			      min_recHitPosition = gemRecHit_loc_x->at(iRecHit);
			      etaP = gemRecHit_etaPartition->at(iRecHit);
			      iEtaRechit = std::to_string(etaP);
			      mupt = mu_propagated_pt->at(iMu);
			      mueta = TMath::Abs(mu_propagated_eta->at(iMu));
			      muphi = mu_propagated_phi->at(iMu);

			      matching = true;
			    }
			  else continue;
			}//end loop rechit

		      if(matching == true)
			{
			  
			  //m_plots["ResidualPlotLocalX_endcappos_dX_"+string]->Fill(min_residual_x);
			  //m_plots["MuPtMatchedLocal_endcappos_dX_"+string+"_iEta_"+iEta]->Fill(mupt);
			  //m_plots["MuEtaMatchedLocal_endcappos_dX_"+string+"_iEta_"+iEta]->Fill(mueta);
			  //m_plots["MuPhiMatchedLocal_endcappos_dX_"+string+"_iEta_"+iEta]->Fill(muphi);

			}
		      else continue;
		    }//end scanning iEta

		} //if positive endcap

	    }//end loop propagated
  
    
		  
		
  }

  
  endJob();
  
}

void gemAnalysis_genParticle::book()
{
  
  m_plots["MuPtPropagated_endcapneg"] = new TH1F("MuPtPropagated_endcapneg",
                                                 "MuPtPropagated; p_{t}; entries",
                                                 66, 0., 100.);

  m_plots["MuPtPropagated_endcappos"] = new TH1F("MuPtPropagated_endcappos",
                                                 "MuPtPropagated; p_{t}; entries",
                                                 66, 0., 100.);

  m_plots["MuEtaPropagated_endcapneg"] = new TH1F("MuEtaPropagated_endcapneg",
                                                  "MuEtaPropagated; #eta; entries",
                                                  60, 0., 3.);

  m_plots["MuEtaPropagated_endcappos"] = new TH1F("MuEtaPropagated_endcappos",
                                                  "MuEtaPropagated; #eta; entries",
                                                  60, 0., 3.);
  
  m_plots["MuPhiPropagated_endcapneg"] = new TH1F("MuPhiPropagated_endcapneg",
                                                  "MuPhiPropagated; #phi; entries",
                                                  120, -3., 3.);
  
  m_plots["MuPhiPropagated_endcappos"] = new TH1F("MuPhiPropagated_endcappos",
                                                  "MuPhiPropagated; #phi; entries",
                                                  120, -3., 3.);

  for(int i=0; i<10; i++)
    {

      float index = 2 - (float)(i*2)/10;

      std::string string = std::to_string(index);

      char buf[16];
      sprintf(buf,"%.1f",index);
      const char* ind = buf;

      char buf1[100] = "ResidualPlotLocalX_endcapneg_dX_" ;
      strcat(buf1,ind);
      const char* name = buf1;

      char buf2[100] = "ResidualPlotLocalX_endcapneg_dX_" ;
      strcat(buf2,ind);
      strcat(buf2,"_endcapneg; residual X ; entries");
      const char* title = buf2;

      m_plots["ResidualPlotLocalX_endcapneg_dX_"+string] = new TH1F(name,
                                                                    title,
                                                                    100, -2., 2.);
    }

  for(int j=1; j<9; j++)
    {
      std::string ieta = std::to_string(j);

      char buffer[16];
      sprintf(buffer,"%d",j);
      const char* j_index = buffer;

      char buffer1[100] = "MuPtPropagated_endcapneg_iEta_";
      strcat(buffer1,j_index);
      const char* name1 = buffer1;
      char buffer3[100] = "MuPtPropagated_endcapneg_iEta_";
      strcat(buffer3,"; p_{t} ; entries");
      const char* title1 = buffer3;
      char buffer2[100] = "MuEtaPropagated_endcapneg_iEta_";
      strcat(buffer2,j_index);
      const char* name2 = buffer2;
      char buffer4[100] = "MuEtaPropagated_endcapneg_iEta_";
      strcat(buffer4,"; #eta ; entries");
      const char* title2 = buffer4;

      m_plots["MuPtPropagated_endcapneg_iEta_"+ieta] = new TH1F(name1,
                                                                title1,
                                                                66,0.,100.);

      m_plots["MuEtaPropagated_endcapneg_iEta_"+ieta] = new TH1F(name2,
                                                                 title2,
                                                                 60,0.,3.);

      for(int i=0; i<10; i++)
        {

          float index = 2 - (float)(i*2)/10;

	  std::string string = std::to_string(index);

          char buf[16];
          sprintf(buf,"%.1f",index);
          const char* ind = buf;

          char buf1[100] = "MuPtMatchedLocal_endcapneg_dX_" ;
          strcat(buf1,ind);
          strcat(buf1,"_iEta_");
          char buf2[100];
          sprintf(buf2,"%d",j);
          const char* iEta = buf2;
          strcat(buf1,iEta);
          const char* name = buf1;

          char buf3[100] = "MuPtMatchedLocal_endcapneg_dX_" ;
          strcat(buf3,ind);
          strcat(buf3,buf2);
          strcat(buf3,iEta);
          strcat(buf3,"; p_{t} ; entries");
          const char* title = buf3;

          m_plots["MuPtMatchedLocal_endcapneg_dX_"+string+"_iEta_"+ieta] = new TH1F(name,
										    title,
										    66, 0., 100.);

        }
    }


}

void gemAnalysis_genParticle::endJob()
{

  m_outFile.cd();
  m_outFile.Write();
  m_outFile.Close();

}

float gemAnalysis_genParticle::round(float var)
{

  float value = (int)(var * 10 + 0.5);
  return (float)value/10;

}

float gemAnalysis_genParticle::dR(float eta1, float eta2, float phi1, float phi2)
{
  auto dp = std::abs(phi1 - phi2);
  auto deta = std::abs(eta1 - eta2);
  if(dp > float(M_PI))
    dp -= float(2 * M_PI);
  float n = TMath::Sqrt(dp*dp + deta*deta);
  return n;
}
