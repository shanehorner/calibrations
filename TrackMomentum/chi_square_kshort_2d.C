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
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

void chi_square_kshort_2d(){

  char name[500];

  sprintf(name, "tgr");

  std::cout << "name check: " << name << std::endl;
  
  gStyle->SetOptStat(0);
  int ngrid = 80;
  int ngridsq = ngrid*ngrid;
  
  float a0start = 0.0178;
  float a0end = 0.0189; 
  float a1start = 0.0;
  float a1end = 0.012;
  
  float da0 = (a0end-a0start)/ngrid;
  float da1 = (a1end-a1start)/ngrid;

  float dbinx = da0/2;
  float dbiny = da1/2;

  float xlo = a0start-dbinx;
  float xhi = a0end+dbinx;

  float ylo = a1start-dbiny;
  float yhi = a1end+dbiny;

  float zmax = 100;

  TH2D *ab = new TH2D("ab", "a vs b weighted X sq.", ngrid,xlo,xhi,ngrid,ylo,yhi);
  TH3D *abx = new TH3D("abx", "(a,b)->X sq.", ngrid,xlo,xhi,ngrid,ylo,yhi,500,0,zmax);
  TGraph2D *g = new TGraph2D();
  
  char findata[500];
  sprintf(findata,"data_tgr/kshort_53800_data.root");
  std::cout << "data input file: " << findata << std::endl;
    
  TFile *fid = TFile::Open(findata);
  
  TGraphErrors *Graph;
  fid->GetObject(name, Graph);

  std::cout << "Graph pointer name: " << Graph->GetName() << " Title: " << Graph->GetTitle() << std::endl;

  int errbinsx = Graph->GetN();
  float dof = (float)errbinsx - 2; //minus 2 parameters: a0, a1.
  
  double wdata[errbinsx];
  double dataerr[errbinsx];
  
  for (int i = 0; i<errbinsx; ++i){
    
    wdata[i] = Graph->GetPointY(i);
    dataerr[i] = Graph->GetErrorY(i);
 
    }
  
  fid->Close();
  
  double chiSquare[ngridsq];
  double chiSquareKeep = 100000;
  double akeep = 100000;
  double a0keep = 100000;
  double bkeep = 100000;
  double a1keep = 100000;
  double ikeep = 100000;
  int problemchild = 0;
  
  for (int i=0; i<ngridsq; ++i){
    bool problem = false;
    int a = i%ngrid;
    int b = i/ngrid;

    float a0 = a0start + (a*da0);
    float a1 = a1start + (b*da1);

    char finsims[500];
    sprintf(finsims, "sim_tgr/kshort_%i.root", i);
    std::cout << "sims input file: " << finsims << std::endl;
    
    TFile *fis = new TFile(finsims);
    TGraphErrors *Graph;
    fis->GetObject(name, Graph);

    double wsims[errbinsx];
    double simserr[errbinsx];
    
    chiSquare[i] = 0;

    std::cout << "a = " << a << " b = " << b << std::endl;

    for (int j=0; j<errbinsx; ++j){
      
      wsims[j] = Graph->GetPointY(j);

      // std::cout << "  wsims[j] = " << wsims[j] << std::endl;
      // std::cout << "  wdata[j] = " << wdata[j] << std::endl;
      
      simserr[j] = Graph->GetErrorY(j);
      if(simserr[j]>5){
	std::cout << " j = " << j << " is a problem child!" << std::endl;
	simserr[j]=5;
	std::cout << "WARNING! decreased simserr[j] to 5." << std::endl;
	problem = true;
      }
      
      double residsq = pow(wdata[j]-wsims[j],2);
      // std::cout << "  resid sq. = " << residsq << std::endl;
      double sigsq = pow(dataerr[j],2)+pow(simserr[j],2);
      // std::cout << "  sigma sq. = " << sigsq << std::endl;
      chiSquare[i]+= (residsq/sigsq);
      // std::cout << "  chisquare update sum = " << chiSquare[i] << std::endl;
      
    }

    if(problem == true){problemchild++;}

    if(std::isinf(chiSquare[i]) || std::isnan(chiSquare[i]) ){
      chiSquare[i]=chiSquare[i-1];
      std::cout << " WARNING !!! INFINITE NAN ENTRY. ERROR. SET TO chiSquare[i-1]: " << std::endl;
    }
      
    if (chiSquare[i]<chiSquareKeep){

      chiSquareKeep = chiSquare[i];
      ikeep = i;
      akeep = a;
      bkeep = b;
      a0keep = a0;
      a1keep = a1;

    }
    
    std::cout << "chisquare entry a: " << a << " b: " << b << " = " << chiSquare[i] << std::endl;
    std::cout << "a0: " << a0 << " , a1: " << a1 << endl;

    g->SetPoint(i,a0,a1,chiSquare[i]);
     
    if(chiSquare[i]<zmax){

      ab->Fill(a0,a1,chiSquare[i]);
      abx->Fill(a0,a1,chiSquare[i]);
      
    }

      
    
    fis->Close();

  }
  TFile *fout = new TFile("contour_plots/ks_ana495_1824_test.root", "recreate");
  TCanvas *cx = new TCanvas("cx", "cx", 100,100,1000,1000);
  g->Draw("surf1");

  std::cout << "minimum raw chi square value:     " << chiSquareKeep << std::endl;
  std::cout << "minimum reduced chi square value: " << (chiSquareKeep/dof) << std::endl;
  std::cout << "respective file number:           " << ikeep << std::endl;
  std::cout << "minimum a: " << akeep << std::endl;
  std::cout << "  so a0 = " << a0keep << std::endl;
  std::cout << "minimum b: " << bkeep << std::endl;
  std::cout << "  so a1 = " << a1keep << std::endl;
  std::cout << "total problem children: " << problemchild << "! Thanks, RooFit." << std::endl;


  double contours[3];
  contours[0] = chiSquareKeep + 1.0;
  contours[1] = chiSquareKeep + 2.30;
  contours[2] = chiSquareKeep + 6.63;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 200,200,900,600);
  c1->Divide(2,1);
  c1->cd(1);
  ab->Draw("box");
  c1->cd(2);
  abx->Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 400,400,1000,1000);
  g->SetName("g");
  g->SetNpx(ngrid);
  g->SetNpy(ngrid);
  g->GetXaxis()->SetTitle("a0 [MS] (GeV)");
  g->GetYaxis()->SetTitle("a1 [DR] (GeV)");
  g->GetZaxis()->SetTitle("#Chi^{2}, 10 d.o.f.");
  g->Write();
  
   //contour plot:
 
  g->GetHistogram()->SetContour(3, contours);
  g->Draw("CONT LIST");
  gPad->Update(); 

  // Grab the contours ROOT created
  TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  if (!conts) { 
    std::cout << "No contours found!" << std::endl;
    return; 
  }

  // Make a clean 2D frame to draw on

  // gROOT->cd();
  
  TCanvas *c = new TCanvas("c","Contours only",800,600);
  TH2F *frame = new TH2F("frame","Contours;a0 [MS] (GeV);a1 [DR] (GeV)",1,0.0178,0.0189,1,0.0,0.012);
  frame->Draw();

  TLatex l;
  l.SetTextSize(0.03);
  char label[32];

  // Loop over contour levels
  std::cout << "conts get size: " << conts->GetSize() << std::endl;
  for (int i=0; i<conts->GetSize(); i++) {
    TList *lcont = (TList*) conts->At(i);
    std::cout << "  lcont get size: " << lcont->GetSize() << std::endl;
    TGraph *gr;
    char name[500];
    sprintf(name, "ctr%i", i);
    for (int j=0; j<lcont->GetSize(); j++) {
      gr = (TGraph*) lcont->At(j);
      gr->SetLineColor(i+1);
      gr->Draw("C");   // draw contour
      char name[500];      // Label contour (take first point as anchor)
      double x0,y0;
      gr->GetPoint(0,x0,y0);
      sprintf(label,"#Delta#Chi^{2}=%.2f", g->GetHistogram()->GetContourLevel(i)-chiSquareKeep);
      l.DrawLatex(x0,y0,label);
      c->Update();
    }
    gr->SetName(name);
    gr->Write();
  }
  // Draw min point
  TMarker *minMarker = new TMarker(a0keep, a1keep, 20);
  minMarker->SetMarkerColor(kBlue);
  minMarker->SetMarkerSize(1.2);
  minMarker->Draw("P");
  minMarker->Write();
  
  // add min label
  TLatex lat;
  char minlabel[500];
  sprintf(minlabel, "  #Chi^{2}min. = %.1f  (%.4f, %.5f)", chiSquareKeep, a0keep, a1keep);
  lat.SetTextSize(0.03);
  lat.DrawLatex(a0keep, a1keep, minlabel);
  
  c->Update();
  c->Write();

  gSystem->Sleep(10000);
  
}
