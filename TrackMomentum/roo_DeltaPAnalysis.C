#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TChain.h"
#include <TTree.h>
#include "TStyle.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TF1.h"
#include <iostream>
#include <map>
#include <set>
#include <TLorentzVector.h>
#include <cmath>
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "RooRandom.h"


void perturbMom(float &p_x, float &p_y, float&p_z, float a0, float a1){

  float pt = sqrt(p_x*p_x + p_y*p_y);
  float dpt = sqrt(pow(a0*pt,2)+pow(a1*pt*pt,2)); //add in quadratures.
  
  TF1* f1 = new TF1 ("f1", "gaus");
  f1->SetParameter(0,1); //scale
  f1->SetParameter(1,0); //mean
  f1->SetParameter(2, dpt); //sigma
  
  // gaussian sample adds a small "dp" perturbing the momentum in each transverse direction.
  p_x = p_x + f1->GetRandom()*sqrt(2); 
  p_y = p_y + f1->GetRandom()*sqrt(2);
  p_z = p_z;
  
}

void roo_DeltaPAnalysis(const int num = 399)
 
{

  // a0 = multiple scattering term, a1 = detector resolution term.
  int ngrid = 80;
  int a = num%ngrid;
  int b = num/ngrid;

  std::cout << "a: " << a << " , b: " << b << endl;

  float a0start = 0.0178;
  float a0end = 0.0189; 
  float a1start = 0.0;
  float a1end = 0.012;

  float da0 = (a0end-a0start)/ngrid;
  float da1 = (a1end-a1start)/ngrid;
    
  float a0 = a0start + (a*da0);
  float a1 = a1start + (b*da1);

  std::cout << "a0: " << a0 << " , a1: " << a1 << endl;

  char outname[500];
  sprintf(outname, "sim_rds/kshort_%d.root", num);
  
  TFile *fout = new TFile(outname, "recreate");

  RooRealVar rPt("rPt", "perturbed pT", 0.0, 5.0);
  RooRealVar rMass("rMass", "perturbed mass", 0.2, 0.7);
  RooArgSet vars(rPt, rMass);
  
  RooDataSet *rmpt = new RooDataSet("rmpt", "m v pt unbin", vars);
  
  for(int i=0;i<10000;++i)
    {
      TChain *ntp = new TChain("ntp_reco_info","ntp");  
     
#include "kshort_variables.C"

      char name[500];
      sprintf(name,"/sphenix/tg/tg01/hf/shanehorner/sims/jet_100m_pythiatest/process_%d.root_kshort.root", i);
     
      std::cout << "Add file " << name << std::endl;
      ntp->Add(name);
      
      int entries = ntp->GetEntries();
     
      std::cout << "entries: " << entries << std::endl;

      for (int i=0; i<entries; ++i){
	ntp->GetEntry(i);
	//if ( crossing1 != crossing2 ){continue;}; //crosscut
	if ( (sqrt(px1*px1+py1*py1) < 0.20) || (sqrt(px2*px2+py2*py2)) < 0.20 ) continue; //ptcut
	if (cosThetaReco < 0.95) continue;
	if ((pow(projected_pathlength_x,2)+pow(projected_pathlength_y,2)) < 0.2) continue;
	if (fabs(projected_pair_dca) > 0.05) continue; //paircut
       	if ((dca3dxy1<0.05) || (dca3dxy2<0.05)) continue;

	// std::cout << "cuts complete. " << endl;

	float decaymass = 0.13957;
	float p_x1 = projected_mom1_x;
	float p_y1 = projected_mom1_y;
	float p_z1 = projected_mom1_z;
	float p_x2 = projected_mom2_x;
	float p_y2 = projected_mom2_y;
	float p_z2 = projected_mom2_z;
	float KSRInvMass = invariant_mass;
	float KSRInvPt = invariant_pt;
	
	//	std::cout << "before: " << "px2: " << p_x2 << " py2: " << p_y2 << " pz2: " << p_z2 << endl;
	//	std::cout << " old invariant mass: " << KSRInvMass << endl;

	if(!(a0 == 0.0 && a1 == 0.0)){

	  perturbMom(p_x1,p_y1,p_z1,a0,a1);
	  perturbMom(p_x2,p_y2,p_z2,a0,a1);
	  
	  //  deltaP(p_x1, p_y1, p_z1, p_x2, p_y2, p_z2, a0, a1);
	  // std::cout << "perturbed succesfully. " << std::endl;

	  
	}


      	// std::cout << "after: " << "px2: " << p_x2 << " py2: " << p_y2 << " pz2: " << p_z2 << endl;
	
	float NewInvMass = 0;
	float NewInvPt = 0;
	double E1 = sqrt(p_x1*p_x1 + p_y1*p_y1 + p_z1*p_z1 + decaymass*decaymass);
	double E2 = sqrt(p_x2*p_x2 + p_y2*p_y2 + p_z2*p_z2 + decaymass*decaymass);
	TLorentzVector v1(p_x1, p_y1, p_z1, E1);
	TLorentzVector v2(p_x2, p_y2, p_z2, E2);
	TLorentzVector tsum;
	tsum = v1 + v2;
	
        NewInvMass = tsum.M();
        NewInvPt = tsum.Pt();

	//     	std::cout << "new inv. mass: " << NewInvMass << " old inv. mass: " << KSRInvMass << endl;
	//	std::cout << "new pt: " << NewInvPt << " old pt: " << KSRInvPt << endl;

	if (NewInvMass<0.2 || NewInvMass>0.7) continue; // perturbed mass cut

	rPt.setVal(NewInvPt);
	rMass.setVal(NewInvMass);
	rmpt->add(vars);
	

	//      std::cout << "new invariant mass: " << NewInvMass << endl;

	
      }
  delete ntp;
    }

  rmpt->Write();
  fout->Close();
  delete fout;
  std::cout << "All Done. " << std::endl;
  
}
