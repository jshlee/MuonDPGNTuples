#define Zmumu_macro_cxx
#include "Zmumu_macro.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TCut.h>
#include <TTree.h>
#include <TString.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLorentzVector.h>
#include <TPad.h>
#include <iostream>
#include <string.h>

#include "vector"

Zmumu_macro::Zmumu_macro(const TString & inFileName,
                         const TString & outFileName) :
  m_inFile(inFileName,"READ"),m_outFile(outFileName,"RECREATE"),fChain(0)
{

  fChain = static_cast<TTree*>(m_inFile.Get("muNtupleProducer/MuDPGTree"));
  Init(fChain);

}

void Zmumu_macro::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   book();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = fChain->LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


     if(mu_nMuons != 0)
       {
	 
	 bool find_mumu_pair = false;
	 //bool muTag_selected_inGE11 = false;
	 bool muProbe_selected_inGE11 = false;
	 //bool muTag_selected_isGEM = false;
	 bool muProbe_selected_isGEM = false;
	 float max_pt_pair = 0, mass_pair = 0;
	 //TLorentzVector muTag_selected;
	 TLorentzVector muProbe_selected;
	 	 
	 for(std::size_t iMuTag = 0; iMuTag < mu_nMuons; iMuTag++)
	   {
	     find_mumu_pair = false;
	     //muTag_selected_inGE11 = false;
	     muProbe_selected_inGE11 = false;
	     //muTag_selected_isGEM = false;
	     muProbe_selected_isGEM = false;
	     mass_pair = 0;
	  	     
	     for(std::size_t iMuProbe = 0; iMuProbe < mu_nMuons; iMuProbe++)
	       {
		 if(iMuTag != iMuProbe && (mu_charge->at(iMuTag))*mu_charge->at(iMuProbe) == -1 && mu_isTight->at(iMuTag) && mu_pt->at(iMuTag) >= 2  && mu_pt->at(iMuProbe) >= 2)
		  		  		   
		   {

		     TLorentzVector muTag;
		     muTag.SetPtEtaPhiM(mu_pt->at(iMuTag),mu_eta->at(iMuTag),mu_phi->at(iMuTag),0.106);

		     TLorentzVector muProbe;
		     muProbe.SetPtEtaPhiM(mu_pt->at(iMuProbe),mu_eta->at(iMuProbe),mu_phi->at(iMuProbe),0.106);
		     
		     TLorentzVector mu_pair = muTag + muProbe;
		     float pt_pair = mu_pair.Pt();

		     if(mu_pair.M() >= 80 && mu_pair.M() <= 100 && pt_pair > max_pt_pair )
		       
		       {
			 max_pt_pair = pt_pair;
			 //mu1_selected.SetPtEtaPhiM(mu_pt->at(iMu1),mu_eta->at(iMu1),mu_phi->at(iMu1),0.106);
			 muProbe_selected.SetPtEtaPhiM(mu_pt->at(iMuProbe),mu_eta->at(iMuProbe),mu_phi->at(iMuProbe),0.106);
			 //mu1_selected_isGEM = mu_isGEM->at(iMu1);
			 muProbe_selected_isGEM = mu_isGEM->at(iMuProbe);
			 mass_pair = mu_pair.M();
			 
			 find_mumu_pair = true;
		       } //if Z mass range
		     else continue;

		   } //if opposite charge
	
	       }//loop for mumu pair
	     
	     //mu1_selected_inGE11 = (TMath::Abs(mu1_selected.Eta()) >= 1.55 && TMath::Abs(mu1_selected.Eta()) <= 2.18);
	     muProbe_selected_inGE11 = (TMath::Abs(muProbe_selected.Eta()) >= 1.55 && TMath::Abs(muProbe_selected.Eta()) <= 2.18);
	     
	     if(find_mumu_pair && muProbe_selected_inGE11)//(mu1_selected_inGE11 || mu2_selected_inGE11) )
	       
	       {
		
		 /*if(mu1_selected_inGE11)
		   {
		     m_plots["MuPt_denominator"]->Fill(mu1_selected.Pt());
		     m_plots["MuEta_denominator"]->Fill(TMath::Abs(mu1_selected.Eta()));
		     m_plots["MuPhi_denominator"]->Fill(mu1_selected.Phi());

		     //for(std::size_t iRecHit = 0; iRecHit < gemRecHit_nRecHits; iRecHit++)
		     //{
		       
		     if(mu1_selected_isGEM)
		       {
			 m_plots["MuPt_numerator"]->Fill(mu1_selected.Pt());
			 m_plots["MuEta_numerator"]->Fill(TMath::Abs(mu1_selected.Eta()));
			 m_plots["MuPhi_numerator"]->Fill(mu1_selected.Phi());
		       }//if matching muon-rechit
		     else continue;  
		     
			 //}//end loop on gem rechit
		     
		   }//if m1 in GE11
		 */
		 //else if(mu2_selected_inGE11)
		 //{
               
		 m_plots["MuPt_denominator"]->Fill(muProbe_selected.Pt());
		 m_plots["MuEta_denominator"]->Fill(TMath::Abs(muProbe_selected.Eta()));
		 m_plots["MuPhi_denominator"]->Fill(muProbe_selected.Phi());
		 
		     //for(std::size_t iRecHit = 0; iRecHit < gemRecHit_nRecHits; iRecHit++)
		       //{
		       
			 if(muProbe_selected_isGEM)
			   {
			     m_plots["MuPt_numerator"]->Fill(muProbe_selected.Pt());
                             m_plots["MuEta_numerator"]->Fill(TMath::Abs(muProbe_selected.Eta()));
                             m_plots["MuPhi_numerator"]->Fill(muProbe_selected.Phi());
			     m_plots["ZMass_success"]->Fill(mass_pair);
			   }//if matching muon-rechit
			 else 
			   {
			     m_plots["ZMass_fail"]->Fill(mass_pair);
			     continue;
			   }
			 //}//end loop on gem rechit
		     
			 //}//if m2 in GE11
		 
			 
	       }//if GEM eta acceptance 
	       else continue;

	   }//end of muon loop

       }//if nMuons != 0
     else continue;
     
   }

   
   m_plots["EfficiencyPt"]->Sumw2();
   m_plots["EfficiencyEta"]->Sumw2();
   m_plots["EfficiencyPhi"]->Sumw2();

   m_plots["EfficiencyPt"]->Divide(m_plots["MuPt_numerator"],m_plots["MuPt_denominator"]);
   m_plots["EfficiencyEta"]->Divide(m_plots["MuEta_numerator"],m_plots["MuEta_denominator"]);
   m_plots["EfficiencyPhi"]->Divide(m_plots["MuPhi_numerator"],m_plots["MuPhi_denominator"]);
   
   //////////////////////////// PLOTS

   TCanvas *c1 = new TCanvas("c1","c1",700,700);
   c1->Divide(1,2);
   
   c1->cd(1);
   gPad->SetGrid();
   gPad->SetPad(0.,0.3,1.,1.);
   gStyle->SetOptStat(0);
      
   m_plots["MuPt_denominator"]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
   m_plots["MuPt_denominator"]->GetXaxis()->SetTitleSize(0.05);
   m_plots["MuPt_denominator"]->GetYaxis()->SetTitle("Entries/1.5 GeV/c bin");
   m_plots["MuPt_denominator"]->GetYaxis()->SetTitleSize(0.05);
   m_plots["MuPt_denominator"]->GetYaxis()->SetTitleOffset(1);
   m_plots["MuPt_denominator"]->SetFillColor(kRed);
   m_plots["MuPt_denominator"]->SetLineColor(kRed);
   m_plots["MuPt_denominator"]->SetFillStyle(3004);
   m_plots["MuPt_denominator"]->Draw("same");
   
   TLatex *t = new TLatex();
   t->SetNDC();
   t->SetTextAlign(33);
   t->SetTextFont(63);
   t->SetTextSizePixels(20);
   t->DrawLatex(.06,.9,"\textbf{CMS} Simulation");

   TPaveText *pt = new TPaveText(.05,.8,.1,.9);
   pt->AddText("Tag & Probe");
   pt->Draw("same");

   gPad->Modified();
   gPad->Update();
   
   m_plots["MuPt_numerator"]->SetFillColor(kRed+2);
   m_plots["MuPt_numerator"]->SetLineColor(kRed+2);
   m_plots["MuPt_numerator"]->SetFillStyle(3005);
   m_plots["MuPt_numerator"]->Draw("same");
   gPad->Modified();
   gPad->Update();

   c1->cd(2);
   gPad->SetGrid();
   gPad->SetPad(0.,0.,1.,0.3);
   m_plots["EfficiencyPt"]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
   m_plots["EfficiencyPt"]->GetXaxis()->SetTitleSize(0.05);
   m_plots["EfficiencyPt"]->GetYaxis()->SetTitle("ratio");
   m_plots["EfficiencyPt"]->GetYaxis()->SetTitleSize(0.05);
   m_plots["EfficiencyPt"]->GetYaxis()->SetTitleOffset(1);
   m_plots["EfficiencyPt"]->Draw();
   gPad->Modified();
   gPad->Update();

   c1->SaveAs("tagprobe.png");
   c1->Clear();
   c1->Close();

   TCanvas *c2 = new TCanvas("c2","c2",700,700);
   
   c2->cd();
   gPad->SetGrid();
   gStyle->SetOptStat(0);

   m_plots["ZMass_success"]->GetXaxis()->SetTitle("M(#mu^{+} #mu^{-}) [GeV/c^{2}]");
   m_plots["ZMass_success"]->GetXaxis()->SetTitleSize(0.04);
   m_plots["ZMass_success"]->GetYaxis()->SetTitle("Entries/0.1 GeV/c^{2} bin");
   m_plots["ZMass_success"]->GetYaxis()->SetTitleSize(0.04);
   m_plots["ZMass_success"]->GetYaxis()->SetTitleOffset(1);
   m_plots["ZMass_success"]->SetFillColor(kRed);
   m_plots["ZMass_success"]->SetLineColor(kBlack);
   m_plots["ZMass_success"]->Draw();

   /*TLatex *t = new TLatex();
   t->SetNDC();
   t->SetTextAlign(33);
   t->SetTextFont(63);
   t->SetTextSizePixels(20);
   t->DrawLatex(.06,.9,"\textbf{CMS} Simulation");*/

   TPaveText *pt2 = new TPaveText(.05,.8,.1,.9);
   TText *text = pt2->AddText("Z #rightarrow #mu^{+}#mu^{-}");
   pt2->Draw("same");

   gPad->Modified();
   gPad->Update();

   c2->SaveAs("invariantmass_success.png");
   c2->Clear();
   c2->Close();
   
   TCanvas *c3 = new TCanvas("c3","c3",700,700);

   c3->cd();
   gPad->SetGrid();
   gStyle->SetOptStat(0);

   m_plots["ZMass_fail"]->GetXaxis()->SetTitle("M(#mu^{+} #mu^{-}) [GeV/c^{2}]");
   m_plots["ZMass_fail"]->GetXaxis()->SetTitleSize(0.04);
   m_plots["ZMass_fail"]->GetYaxis()->SetTitle("Entries/0.1 GeV/c^{2} bin");
   m_plots["ZMass_fail"]->GetYaxis()->SetTitleSize(0.04);
   m_plots["ZMass_fail"]->GetYaxis()->SetTitleOffset(1);
   m_plots["ZMass_fail"]->SetFillColor(kRed);
   m_plots["ZMass_fail"]->SetLineColor(kBlack);
   m_plots["ZMass_fail"]->Draw();

   /*TLatex *t = new TLatex();
   t->SetNDC();
   t->SetTextAlign(33);
   t->SetTextFont(63);
   t->SetTextSizePixels(20);
   t->DrawLatex(.06,.9,"\textbf{CMS} Simulation");*/

   TPaveText *pt3 = new TPaveText(.05,.8,.1,.9);
   TText *text2 = pt3->AddText("Z #rightarrow #mu^{+}#mu^{-}");
   pt3->Draw("same");

   gPad->Modified();
   gPad->Update();

   c3->SaveAs("invariantmass_fail.png");
   c3->Clear();
   c3->Close();

   endJob();
   
}

void Zmumu_macro::book()
{

  m_plots["ZMass_fail"] = new TH1F("ZMass_fail",
			      "",
			      200,80.,100.);
  
  m_plots["ZMass_success"] = new TH1F("ZMass_success",
				      "",
				      200,80.,100.);


  m_plots["MuPt_denominator"] = new TH1F("MuPt_denominator",
					 "; p_{t}; entries",
					 66,0.,100.);

  m_plots["MuEta_denominator"] = new TH1F("MuEta_denominator",
					  "; #eta; entries",
					  60, 0.,3.);

  m_plots["MuPhi_denominator"] = new TH1F("MuPhi_denominator",
					  "; #phi; entries",
                                          70, -3.5,3.5);
    
  m_plots["MuPt_numerator"] = new TH1F("MuPt_numerator",
				       "; p_{t}; entries",
				       66,0.,100.);

  m_plots["MuEta_numerator"] = new TH1F("MuEta_numerator",
					"; #eta; entries",
					60, 0.,3.);

  m_plots["MuPhi_numerator"] = new TH1F("MuPhi_numerator",
					"; #phi; entries",
					70, -3.5,3.5);

  m_plots["EfficiencyPt"] = new TH1F("EfficiencyPt",
				     "; p_{t}; efficiency",
				     66, 0.,100.);
  
  m_plots["EfficiencyEta"] = new TH1F("EfficiencyEta",
				      "; #eta; entries",
				      60, 0.,3.);
  
  m_plots["EfficiencyPhi"] = new TH1F("EfficiencyPhi",
				      "; #phi; entries",
				      70, -3.5,3.5);

}


void Zmumu_macro::endJob()
{

  m_outFile.cd();
  m_outFile.Write();
  m_outFile.Close();

}
