#include <TFile.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TF1.h>
#include <TChain.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraphPainter.h>
#include <TPaveStats.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>



void compare_tgr_kshort_simvdata(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  // input file number:
  
  int num = 3634;

  char outname[500];
  sprintf(outname, "tgr_data_sim_compare/ana_ks_53800_1824_0.2_1.8_3.5c2.root");

 

  //switches: roo vs root

  char name[500];

  sprintf(name, "tgr"); 
  
  TCanvas *c1 = new TCanvas("c1", "", 200,200,1200,1200);
  char simfile[500];
  sprintf(simfile, "sim_tgr/kshort_%i.root", num);
  TFile *fis = new TFile(simfile);

  TGraphErrors *grs = (TGraphErrors*) fis->Get(name);
  if(!grs){std::cout << "failed to get input file fis Graph" << std::endl;
    return;}
  grs->SetDrawOption("AP");
  grs->SetMarkerStyle(20);
  grs->SetMarkerColor(kBlue);
  grs->SetFillStyle(0);
 

  TFile *fid = new TFile("data_tgr/kshort_53800_data.root");

  TGraphErrors *grd = (TGraphErrors*) fid->Get(name);
  if(!grd){std::cout << "failed to get input file fid Graph" << std::endl;
    return;}
  grd->SetDrawOption("AP");
  grd->SetMarkerStyle(20);
  grd->SetMarkerColor(kRed);
  grd->SetMarkerSize(1);
  //grd->Fit(dptd); 
  

  // fid->Close();

  TMultiGraph *mg = new TMultiGraph();
  
  mg->SetTitle("Sim (blue), data (red)");
  mg->SetName("mg");
  mg->Add(grs);
  mg->Add(grd);
  mg->GetXaxis()->SetTitle("P_{T} (GeV)");
  mg->GetYaxis()->SetTitle("K-Short mass width (MeV)");
  mg->Draw("AP");

 TFile *fout = new TFile(outname, "recreate");
  
  mg->Write("mg");
  fout->Close();
  


}
