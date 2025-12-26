#include <TFile.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooMsgService.h>

using namespace RooFit;

void single_plot_analyze_kshort_mom()
{
  //Useful : set DOF, pT range, mass range of interest.
  static const int nslices = 12;
  double ptmax = 1.8;
  double ptmin = 0.2;
  double mlo = 0.43;
  double mhi = 0.56;

  //Choose data/sim. For simulation input, set file number:
  int num = 70;
  bool dataSwitch = true;
  bool simSwitch = !dataSwitch;

  gStyle->SetOptStat(0);

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  char finName[500];
  char foutName[500];
  if (simSwitch)
    {
      sprintf(finName, "sim_rds/kshort_%i.root", num);
      sprintf(foutName, "sim_tgr/kshort_test.root");
    }
  if (dataSwitch)
    {
      sprintf(finName, "data_rds/run53877_0_1824_kshort_ana495_ch_cut.root");
      sprintf(foutName, "data_tgr/kshort_53800_data.root");
    }
 
  std::cout << "infile : " << finName << " outfile: " << foutName << std::endl;
      
  TFile *fin = new TFile(finName);
  
  if (dataSwitch)
    {
      
    }
  RooDataSet *rmpt = (RooDataSet*) fin->Get("rmpt");
  RooRealVar *rMass = (RooRealVar*) rmpt->get()->find("rMass"); // mass
  
  RooRealVar mu("mu", "mu", 0.497, 0.47, 0.51); // mean
  RooRealVar sigma("sigma", "sigma", 0.007, 0.001, 0.03); // stdev
  RooRealVar lambda("lambda", "lambda",0,-20,20);  //decay param.
  RooRealVar nsig("nsig", "nsig", 100000, 0, 500000); // no. entries in signal
  RooRealVar nbkg("nbkg", "nbkg", 1000, 0, 100000);  // no. entries background
  RooGaussian gauss("gauss", "gauss", *rMass, mu, sigma);
  RooExponential expo("expo", "expo", *rMass, lambda);

  RooAddPdf model("model", "model", RooArgList(gauss,expo), RooArgList(nsig,nbkg));
      
  double width[nslices];
  double ptval[nslices];
  double werror[nslices];
  double pterr[nslices];

  double dslice = (ptmax - ptmin)/nslices;

  rmpt->Print();

  TCanvas *csl = new TCanvas ("csl", "slices", 50,50,1000,800);
  int cslicey = sqrt(nslices);
  int cslicex = nslices/cslicey;
  csl->Divide(cslicex,cslicey);

  for(int i = 0; i < nslices; ++i)
    {
      double slo = ptmin + dslice * (double) i;
      double shi = ptmin + dslice * (double) (i+1);
      ptval[i] = 0.5 * (slo + shi);
      pterr[i] = 0.0;

      char ptslice[500];
      sprintf(ptslice, "rPt >= %.6f && rPt < %.6f", slo, shi);	  
	  
      RooDataSet* sliceData = (RooDataSet*) rmpt->reduce(Cut(ptslice));
      model.fitTo(*sliceData, PrintLevel(-1));

      double provmean = sliceData->mean(*rMass);
      double provsig = sliceData->sigma(*rMass);

      std::cout << "provmean: " << provmean << std::endl;

      rMass->setRange(mlo,mhi);
      
      csl->cd();
      csl->cd(i+1);
      

      RooPlot *frame = rMass->frame();
      frame->SetXTitle("MeV/(c^2)");
      char title[500];
      sprintf(title,"Mass, pT %.2f - %.2f",slo,shi);
      frame->SetTitle(title);
      
      sliceData->plotOn(frame);
      model.plotOn(frame, LineColor(kBlue));
      frame->Draw();
      csl->Update();

      width[i] = 1000.0 * sigma.getVal();
      werror[i] = 1000.0 * sigma.getError();

      delete sliceData;
	  
    }
  csl->Update();
  TFile *fout = new TFile(foutName,"recreate");

  TCanvas *ctg = new TCanvas("ctg", "pt -> mass width, error", 100,20,600,800);
  ctg->cd();
  gPad->SetLeftMargin(0.15);
  
  TGraphErrors *tgr = new TGraphErrors(nslices, ptval, width, pterr, werror);
  tgr->SetMinimum(0.7*width[0]);
  tgr->SetMaximum(1.3*width[nslices-1]);
  tgr->SetMarkerStyle(20);
  tgr->SetLineColor(kBlue);
  tgr->SetMarkerSize(1);
  tgr->SetName("tgr");
  tgr->SetTitle("RooFit Results");
  tgr->GetXaxis()->SetTitle("pT (GeV/c)");
  tgr->GetYaxis()->SetTitle("Ks mass width (MeV)");
  tgr->Draw("AP");
  tgr->Write();

  fout->Close();
      

  
}
