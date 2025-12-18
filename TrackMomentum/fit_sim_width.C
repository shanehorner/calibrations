#include <TFile.h>
#include <TNtuple.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooMsgService.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace RooFit;

void fit_sim_width()
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  double ptmin = 0.4;
  double ptmax = 2;
  int nslices = 12;
  double rResMin = -0.3;
  double rResMax = 0.3;
  RooRealVar rPt("rPt", "rPt", ptmin,ptmax);
  RooRealVar rRes("rRes", "rRes", rResMin,rResMax);
  RooArgSet vars(rPt, rRes);
  RooDataSet *data = new RooDataSet("data", "data", vars);

  long candidates = 0;
  for(int process = 0; process<50; ++process)
    {
      TChain *ntp = new TChain("ntp_track");
      float pt, gpt, quality, geta, gvz, gtrackID, gprimary, gflavor, gnmaps;
      ntp->SetBranchAddress("geta", &geta);
      ntp->SetBranchAddress("gvz", &gvz);
      ntp->SetBranchAddress("pt", &pt);
      ntp->SetBranchAddress("gpt", &gpt);
      ntp->SetBranchAddress("quality", &quality);
      ntp->SetBranchAddress("gtrackID", &gtrackID);
      ntp->SetBranchAddress("gprimary", &gprimary);
      ntp->SetBranchAddress("gflavor", &gflavor);
      ntp->SetBranchAddress("gnmaps", &gnmaps);

      char name[500];
      sprintf(name, "/sphenix/tg/tg01/hf/shanehorner/sims/jet_100m_pythiatest/process_%i_g4svtx_eval.root", process);

      ntp->Add(name);
      int entries = ntp->GetEntries();
      std::cout << "process " << process << " entries: " << entries << std::endl;
      for(int entry = 0; entry<entries; ++entry)
	{
	  ntp->GetEntry(entry);
	  if (quality > 5) continue;
	  if (gtrackID<0) continue;
	  if (gprimary != 1) continue;
	  if (abs(gvz)>10) continue;
	  if (abs(gflavor) != 211) continue;
	  if (gnmaps < 3) continue;
	  if (abs(geta) > 1) continue;
	  if (gpt < ptmin) continue;
	  if (gpt > ptmax) continue; 
	  float dpt = (pt-gpt)/gpt;
	  if (dpt > rResMax) continue;
	  if (dpt < rResMin) continue;
	  rRes.setVal(dpt);
	  rPt.setVal(gpt);
	  data->add(vars);
	  
	}
      
      delete ntp;
      
    }
  std::cout << "total candidates: " << candidates << std::endl;
  
  RooRealVar mu("mu", "mu", 0, -0.07, 0.07); // mean
  RooRealVar sigma("sigma", "sigma", 0.01, 0.0001, 0.02); // stdev
  RooRealVar lambda("lambda", "lambda",0,-20,20);  //decay param.
  RooRealVar nsig("nsig", "nsig", 100000, 0, 1000000); // no. entries in signal
  RooRealVar nbkg("nbkg", "nbkg", 1000, 0, 100000);  // no. entries background 
  RooRealVar alpha("alpha", "alpha", 1.5, 0.5, 5);
  RooRealVar n("n", "n", 2, 0.1, 10);

  RooGaussian gauss("gauss", "gauss", rRes, mu, sigma);
  RooExponential expo("expo", "expo", rRes, lambda);

  RooAddPdf model("model", "model", RooArgList(gauss,expo), RooArgList(nsig,nbkg));

  RooRealVar muwide("muwide", "muwide", 0, -0.1,0.1);
  RooRealVar sigwide("sigwide", "sigwide", 0.04,0.02,0.4);
  RooGaussian gausswide("gausswide", "gausswide", rRes, muwide, sigwide);

  RooAddPdf doubleGauss("doubleGauss","double Gauss", RooArgList(gauss,gausswide), RooArgList(nsig,nbkg));
   
  
  RooCBShape cb("cb", "Crystal Ball PDF", rRes, mu, sigma, alpha, n);

       
  double width[nslices];
  double ptval[nslices];
  double werror[nslices];
  double pterr[nslices];

  double dslice = (ptmax-ptmin)/nslices;
  
  TCanvas *csl = new TCanvas ("csl", "slices", 50,50,1000,800);
  int cslicey = sqrt(nslices);
  int cslicex = nslices/cslicey;
  csl->Divide(cslicex,cslicey);

  double wmax = 0;
  double wmin = 100;
  
  for(int slice = 0; slice < nslices; ++slice)
    {
      double slo = ptmin + dslice * (double) slice;
      double shi = ptmin + dslice * (double) (slice + 1);
      ptval[slice] = 0.5 * (slo + shi);
      pterr[slice] = 0.0;

      char slicecut[500];
      sprintf(slicecut, "rPt >= %.6f && rPt < %.6f", slo, shi);
      
      RooDataSet* sliceData = (RooDataSet*) data->reduce(Cut(slicecut));
      mu.setVal(sliceData->mean(rRes));
      sigma.setVal(sliceData->sigma(rRes));

      double rlo = mu.getVal()-3*sigma.getVal();
      double rhi = mu.getVal()+3*sigma.getVal();

      rRes.setRange("tight", rlo, rhi);

      gauss.fitTo(*sliceData, Range("tight"), PrintLevel(-1));
      csl->cd();
      csl->cd(slice+1);

      RooPlot *frame = rRes.frame(Range("tight"));
      frame->SetXTitle("(#DeltaP_{T})/P_{T}");
      char title[500];
      sprintf(title, "#DeltaP_{T}/P_{T}, %.2f - %.2f GeV/c",slo,shi);
      frame->SetTitle(title);

      sliceData->plotOn(frame);
      model.plotOn(frame, Range("tight"), NormRange("tight"), LineColor(kBlue));
      //doubleGauss.plotOn(frame, LineColor(kBlue));
      //gauss.plotOn(frame, LineColor(kRed), LineStyle(kDashed));
      //gausswide.plotOn(frame, LineColor(kGreen), LineStyle(kDashed));
      frame->Draw();
      csl->Update();

      width[slice] = 1000.0 * sigma.getVal();
      werror[slice] = 1000.0 * sigma.getError();

      std::cout << "slice " << slice << " dpt/pt: " << std::endl;
      std::cout << width[slice] << " +/- " << werror[slice] << std::endl;

      if (wmax < width[slice]) wmax = 1.1*width[slice];
      if (wmin > width[slice]) wmin = 0.9*width[slice];
      
      delete sliceData;

    }
  csl->Update();
  
  TGraphErrors *tgr = new TGraphErrors(nslices, ptval, width, pterr, werror);
  tgr->SetMinimum(wmin);
  tgr->SetMaximum(wmax);
  tgr->SetMarkerStyle(20);
  tgr->SetLineColor(kBlue);
  tgr->SetMarkerSize(1);
  tgr->SetName("tgr");
  tgr->SetTitle("RooFit Results");
  tgr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  tgr->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T} (MeV)");

  TF1 *f = new TF1("f", "sqrt(([0])^2+([1]*x)^2)", ptmin, ptmax);
  tgr->Fit(f, "R");

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  tgr->Draw("AP");
  f->Draw("same");
  
}
