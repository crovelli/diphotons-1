#include "TString.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>

const TString lumiString = "CMS Preliminary, #sqrt{s}=13 TeV, #intLdt=1.9 fb^{-1}";

const int nEtaBins = 8;   
const double etaBinLimits[nEtaBins+1]  = {-2.5, -2., -1.5, -1., 0., 1., 1.5, 2., 2.5};
const double etaBinCenters[nEtaBins]   = {-2.25, -1.75, -1.25, -0.5, 0.5, 1.25, 1.75, 2.25};
const double etaBinHalfWidth[nEtaBins] = { 0.25,  0.25,  0.25, 0.5, 0.5, 0.25, 0.25, 0.25};

// ----------------------------------------
// Data efficiencies and statistical errors
double data[nEtaBins] = {
  7.89948e-01, 7.49984e-01, 8.62455e-01, 8.59384e-01, 8.62066e-01, 8.65602e-01, 7.61567e-01, 7.95977e-01 
};

// statistical only errors 
double dataErrStat[nEtaBins] = {
  1.56577e-03, 1.41776e-03, 9.96317e-04, 5.77391e-04, 5.66172e-04, 1.00516e-03, 1.51788e-03, 1.50448e-03
};

// MC efficiencies and errors - C&C
double mc[nEtaBins] = {
  0.724601, 0.759899, 0.892286, 0.889741, 0.888193, 0.890842, 0.776299, 0.756487
};

// statistical only errors 
double mcErr[nEtaBins] = {
  0.,0.,0.,0.,0.,0.,0.,0.
};

void drawResultsEta(){

  // Stat error
  double dataErr[nEtaBins];

  std::cout << "only statistical error" << std::endl;
  for (int ii=0; ii<nEtaBins; ii++ ) {
    dataErr[ii] = dataErrStat[ii];
  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << "scale factors" << std::endl;

  double sf[nEtaBins];
  double sfErrTot[nEtaBins];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    sf[iEta] = data[iEta]/mc[iEta];
    float sigmaDoD   = dataErr[iEta]/data[iEta];
    float sigmaMCoMC = mcErr[iEta]/mc[iEta];
    sfErrTot[iEta] = sf[iEta]*sqrt( (sigmaDoD*sigmaDoD) + (sigmaMCoMC*sigmaMCoMC) );
    std::cout << sf[iEta] << " +/- " << sfErrTot[iEta] << std::endl;
  }

  // Draw all canvases
  TString cname = "sfEffEta";
  TCanvas *c1 = new TCanvas(cname, cname, 10,10,700,700);
  c1->SetFillColor(kWhite);
  c1->Draw();
  TPad *pad1 = new TPad("main","",0, 0.3, 1.0, 1.0);
  pad1->SetTopMargin(0.20);
  pad1->SetBottomMargin(0.02);
  pad1->SetGrid();
  TPad *pad2 = new TPad("ratio", "", 0, 0, 1.0, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.30);
  pad2->SetGrid();
  
  pad1->Draw();
  pad2->Draw();
  
  // Create and fill arrays for graphs for this eta bin
  double *dataSlice    = new double[nEtaBins];
  double *dataSliceErr = new double[nEtaBins];
  double *mcSlice      = new double[nEtaBins];
  double *mcSliceErr   = new double[nEtaBins];
  double *sfSlice      = new double[nEtaBins];
  double *sfSliceErr   = new double[nEtaBins];
  for(int ieta = 0; ieta<nEtaBins; ieta++){
    dataSlice   [ieta] = data     [ieta];
    dataSliceErr[ieta] = dataErr  [ieta];
    mcSlice     [ieta] = mc       [ieta];
    mcSliceErr  [ieta] = mcErr    [ieta];
    sfSlice     [ieta] = sf       [ieta];
    sfSliceErr  [ieta] = sfErrTot [ieta];
  }
  
  // Create and configure the graphs   
  TGraphErrors *grData = new TGraphErrors(nEtaBins, etaBinCenters, dataSlice, etaBinHalfWidth, dataSliceErr);
  TGraphErrors *grMc   = new TGraphErrors(nEtaBins, etaBinCenters, mcSlice, etaBinHalfWidth, mcSliceErr);
  TGraphErrors *grSf   = new TGraphErrors(nEtaBins, etaBinCenters, sfSlice, etaBinHalfWidth, sfSliceErr);
  
  grData->SetLineColor(kBlack);
  grData->SetMarkerColor(kBlack);
  grData->SetMarkerStyle(20);
  grData->SetMarkerSize(1.);
  
  int ci = TColor::GetColor("#99ccff");
  grMc->SetFillColor(kGreen-8);
  ci = TColor::GetColor("#3399ff");
  grMc->SetLineColor(kGreen+4);
  grMc->SetMarkerStyle(22);
  grMc->SetMarkerColor(kGreen+4);
  grMc->SetMarkerSize(1.);
  
  ci = TColor::GetColor("#99ccff");
  grSf->SetFillColor(kGreen-8);
  ci = TColor::GetColor("#3399ff");
  grSf->SetLineColor(kGreen+4);
  grSf->SetMarkerStyle(20);
  grSf->SetMarkerColor(kGreen+4);
  grSf->SetMarkerSize(1.);
  
  // Create and configure the dummy histograms on which to draw the graphs
  TH2F *h1 = new TH2F("dummy1","", 100, -2.6, 2.6, 100, 0.6, 1.1);
  h1->GetYaxis()->SetTitle("Efficiency");
  h1->SetStats(0);
  h1->GetXaxis()->SetLabelSize(0);
  h1->GetXaxis()->SetNdivisions(505);
  h1->GetXaxis()->SetDecimals();
  h1->GetYaxis()->SetTitleOffset(0.8);
  h1->GetYaxis()->SetTitleSize(0.05);
  TH2F *h2 = new TH2F("dummy2","", 100, -2.6, 2.6, 100, 0.8, 1.2);
  h2->GetXaxis()->SetTitle("#eta");
  h2->GetYaxis()->SetTitle("Scale Factor");
  h2->GetXaxis()->SetTitleOffset(1.0);
  h2->GetXaxis()->SetTitleSize(0.1);
  h2->GetYaxis()->SetTitleOffset(0.4);
  h2->GetYaxis()->SetTitleSize(0.1);
  h2->GetXaxis()->SetLabelSize(0.08);
  h2->GetYaxis()->SetLabelSize(0.08);
  h2->GetYaxis()->SetNdivisions(505);
  h2->GetYaxis()->SetDecimals();
  h2->SetStats(0);
  
  TLegend *leg = new TLegend(0.65,0.1,0.9,0.25);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(grData, "Data", "pl");
  leg->AddEntry(grMc, "Simulation DY", "pFlE");
  
  TLatex *latLumi = new TLatex(0, 0, lumiString);
  
  // --------------------------------------
  // Draw the efficiencies
  pad1->cd();
  h1->Draw();
  grMc  ->Draw("2same");
  grMc  ->Draw("pZ,same");
  grData->Draw("PEZ,same");
  leg->Draw("same");
  latLumi->Draw("same");
  // Draw the scale factors
  pad2->cd();
  h2->Draw();
  grSf  ->Draw("2same");
  grSf  ->Draw("pEZ,same");
  // Save into a file
  TString fname = cname;
  fname += ".pdf";
  c1->Print(fname);
}



