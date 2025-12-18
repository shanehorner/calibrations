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

void fast_bulk_plot_analyze_kshort_mom(int num = 100)
{
  //Useful : set DOF, pT range, mass range of interest.
  static const int nslices = 12;
  double ptmax = 1.8;
  double ptmin = 0.2;
  double mlo = 0.43;
  double mhi = 0.56;

  gStyle->SetOptStat(0);

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  char finName[500];
  sprintf(finName, "sim_rds/kshort_%i.root", num);
  char foutName[500];
  sprintf(foutName, "sim_tgr/kshort_%i.root", num);

  std::cout << "infile : " << finName << " outfile: " << foutName << std::endl;
      
  TFile *fin = new TFile(finName);

  RooDataSet *rmpt = (RooDataSet*) fin->Get("rmpt");
  RooRealVar *rMass= (RooRealVar*) rmpt->get()->find("rMass"); // mass
  rMass->setRange(mlo,mhi);
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

      width[i] = 1000.0 * sigma.getVal();
      werror[i] = 1000.0 * sigma.getError();

      delete sliceData;
	  
    }

  TFile *fout = new TFile(foutName,"recreate");

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
  tgr->Write();

  fin->Close();
  fout->Close();
  
}
