#include <TFile.h>
#include <TNtuple.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace RooFit;

void multirun_analyze_kshort_mom()
{

  char outname[500];
  int runs[] = {53876, 53877, 53879};
  int segstart = 0;
  int nseg[] = {2483, 2627, 1638};

  sprintf(outname, "data_rds/run53800_all_kshort_ana495.root");
  TFile *fout = new TFile(outname, "recreate");

  double mlo = 0.35;
  double mhi = 0.65;
  double ptlo = 0.2;
  double pthi = 5.0;

  RooRealVar rMass("rMass", "rMass", mlo, mhi);
  RooRealVar rPt("rPt", "rPt", ptlo, pthi);
  RooArgSet vars(rPt,rMass);
  RooDataSet *rmpt = new RooDataSet("rmpt", "roo: (pt, mass)", vars);

  int total = 0;

  for (int runIdx = 0; runIdx < 3; ++runIdx)
    {
      int run = runs[runIdx];
      int nsegs = nseg[runIdx];


      for (int segment = segstart; segment < segstart+nsegs; ++segment)
	{
	  int skip = 0;
	  if (segment < 70) skip = 1000;
	  else skip = 2500;

	  int max = 10000/skip;

	  int candidates = 0;
     
	 
	  TChain *ntp = new TChain("ntp_reco_info");

#include "kshort_variables.C"
      
	  char name[500];
	  sprintf(name, "/sphenix/tg/tg01/hf/frawley/data/ana495_includes_kshort_run53800/_%i-%i.root_kshort.root",run, segment);
	  ntp->Add(name);

	  int entries = ntp->GetEntries();

	  for (int i=0; i<entries; ++i)
	    {
	      ntp->GetEntry(i);
	      if (invariant_pt < ptlo) continue;
	      if ((invariant_mass < mlo) || (invariant_mass > mhi)) continue;
	      //if (crossing1 != crossing2) continue; //crosscut
	      if (sqrt(px1*px1+py1*py1) < 0.20) continue;
	      if (sqrt(px2*px2+py2*py2) < 0.20) continue; //ptcut
	      if (cosThetaReco < 0.95) continue;
	      if ((pow(projected_pathlength_x,2)+pow(projected_pathlength_y,2)) < 0.2) continue;
	      if (fabs(projected_pair_dca) > 0.05) continue; //paircut
	      if (charge1 == charge2) continue;
	      // if( track1_mvtx_states < 3 ) continue;
	      // if( track2_mvtx_states < 3 ) continue; 
	      if( abs(dca3dxy1) < 0.0 ) continue;
	      if( abs(dca3dxy2) < 0.0 ) continue;

	      // std::cout << "cuts complete. " << endl;

	      //std::cout << " m: " << invariant_mass << " pt: " << invariant_pt << std::endl;

	      rMass.setVal(invariant_mass);
	      rPt.setVal(invariant_pt);

	      rmpt->add(vars);
	      ++candidates;

	    }
      
	  total += candidates;
	  std::cout << "run " << run << " segment " << segment <<  " entries: " << entries << std::endl;
	  std::cout << candidates << " kshort candidates" << std::endl;

	  delete ntp;

	}

    }
  std::cout << total << " total kshort candidates. " << std::endl;
  rmpt->Write();
  fout->Close();

  std::cout << "All done. " << std::endl;
  
}
