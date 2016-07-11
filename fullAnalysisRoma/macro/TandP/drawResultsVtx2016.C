#include "TString.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>

const int nEtaBins = 1;

const TString lumiString = "CMS Preliminary, #sqrt{s}=13 TeV, #intLdt=4.0 fb^{-1}";

// EB
const int nNvtxBinsEB = 6;   
const double nvtxBinLimitsEB[nNvtxBinsEB+1]  = { 0., 5., 10., 15., 20., 25., 50. };
const double nvtxBinCentersEB[nNvtxBinsEB]   = { 2.5, 7.5, 12.5, 17.5, 22.5, 37.5 };
const double nvtxBinHalfWidthEB[nNvtxBinsEB] = { 2.5, 2.5, 2.5, 2.5, 2.5, 12.5 };
const TString etaLimitsStringArrayEB[nEtaBins] = { "0. < |#eta| < 1.4442" };

// EE
const int nNvtxBinsEE = 6;   
const double nvtxBinLimitsEE[nNvtxBinsEE+1]  = { 0., 5., 10., 15., 20., 25., 50. };
const double nvtxBinCentersEE[nNvtxBinsEE]   = { 2.5, 7.5, 12.5, 17.5, 22.5, 37.5 };
const double nvtxBinHalfWidthEE[nNvtxBinsEE] = { 2.5, 2.5, 2.5, 2.5, 2.5, 12.5 };
const TString etaLimitsStringArrayEE[nEtaBins] = { "1.566 < |#eta| < 2.5" };


// ----------------------------------------
// Data efficiencies and statistical errors
double dataEB[nEtaBins][nNvtxBinsEB] = {
  { 
    8.54377e-01, 8.55774e-01, 8.59374e-01, 8.64640e-01, 8.68335e-01, 8.70869e-01
  }
};

double dataEE[nEtaBins][nNvtxBinsEE] = {
  { 
    // ALL EE
    // 3.71575e-01, 7.05957e-01, 7.63387e-01, 7.96578e-01, 8.09098e-01, 8.16657e-01

    // EE, eta<2.1
    // 3.11080e-01, 6.88481e-01, 7.46353e-01, 7.85390e-01, 7.99473e-01, 8.17624e-01
    
    // ALL EE, N-1 phiso
    // 9.50878e-01, 9.38365e-01, 9.30607e-01, 9.22462e-01, 9.11179e-01, 8.99571e-01

    // ALL EE, full sel, vs RHO 
    // 6.31302e-01, 7.51349e-01, 7.93608e-01, 8.07686e-01, 8.22531e-01, 8.22953e-01
    
    // ALL EE, corr vs pT only, vs RHO
    8.93515e-01, 8.33850e-01, 7.59603e-01, 6.77719e-01, 6.15368e-01, 5.36136e-01 
  }
};

// statistical only errors 
double dataErrStatEB[nEtaBins][nNvtxBinsEB] = {
  { 
    6.06425e-03, 9.46696e-04, 6.26298e-04, 3.71647e-04, 1.13609e-03, 4.06263e-04 
  }
};

double dataErrStatEE[nEtaBins][nNvtxBinsEE] = {
  { 
    // |eta|<2.5    
    // 1.23516e-02, 2.38062e-03, 1.32544e-03, 1.40598e-03, 2.39583e-03, 5.40049e-03 

    // |eta|<2.1  
    // 3.58185e-02, 3.03384e-03, 1.63594e-03, 1.42734e-03, 2.94480e-03, 6.40698e-03

    // ALL EE, N-1 phiso 
    // 8.32998e-03, 1.85787e-03, 9.96354e-04, 8.35764e-04, 1.41344e-03, 4.49701e-03

    // ALL EE, full sel, vs RHO 
    // 3.07043e-03, 1.17834e-03, 9.68440e-04, 2.02306e-03, 4.79005e-03, 1.27423e-02

    // ALL EE, corr vs pT only, vs RHO
    2.63218e-03, 8.99174e-04, 1.19687e-03, 2.21238e-03, 5.37082e-03, 1.50387e-02
  }
};

// ----------------------------------------
// MC efficiencies and errors - C&C
double mcEB[nEtaBins][nNvtxBinsEB] = {
  { 
    0.879511, 0.880671, 0.885803, 0.892164, 0.898139, 0.904226
  }
};

double mcEE[nEtaBins][nNvtxBinsEE] = {
  { 
    // |eta|<2.5 
    // 0.251313, 0.624328, 0.742932, 0.78909, 0.81591, 0.83126

    // |eta|<2.1
    // 0.219307, 0.624941, 0.750271, 0.795384, 0.823233, 0.841198

    // ALL EE, N-1 phiso
    // 0.956275, 0.952662, 0.946172, 0.939112, 0.931795, 0.922229

    // ALL EE, full sel, vs RHO 
    // 0.560185, 0.744796, 0.803245, 0.827278, 0.836836, 0.836765

    // ALL EE, corr vs pT only, vs RHO
    0.894411, 0.832542, 0.763507, 0.695146, 0.634624, 0.570544
  }
};

// statistical only errors - 218pb
double mcErrEB[nEtaBins][nNvtxBinsEB] = {
  { 
    0.,0.,0.,0.,0.,0.
  }
};

double mcErrEE[nEtaBins][nNvtxBinsEE] = {
  { 
    0.,0.,0.,0.,0.,0.
  }
};


void drawResults(){

  // stat error only
  double dataErrEB[nEtaBins][nNvtxBinsEB];
  double dataErrEE[nEtaBins][nNvtxBinsEE];

  std::cout << "only statistical error" << std::endl;
  for (int ii=0; ii<nNvtxBinsEB; ii++ ) {
    dataErrEB[0][ii] = dataErrStatEB[0][ii];
  }
  for (int ii=0; ii<nNvtxBinsEE; ii++ ) {
    dataErrEE[0][ii] = dataErrStatEE[0][ii];
  }

  std::cout << "================================" << std::endl;
  std::cout << "EB" << std::endl;
  for (int ii=0; ii<nNvtxBinsEB; ii++ ) 
    std::cout << ii << ", nominal = " << dataEB[0][ii] << ", statErr = " << dataErrStatEB[0][ii] << std::endl; 
  std::cout << "================================" << std::endl;
  std::cout << "EE" << std::endl;
  for (int ii=0; ii<nNvtxBinsEE; ii++ ) 
    std::cout << ii << ", nominal = " << dataEE[0][ii] << ", statErr = " << dataErrStatEE[0][ii] << std::endl; 
  std::cout << "================================" << std::endl;


  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << "scale factors: EB" << std::endl;
  // Scale factors and errors
  double sfEB[nEtaBins][nNvtxBinsEB];
  double sfErrTotEB[nEtaBins][nNvtxBinsEB];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iNvtx=0; iNvtx<nNvtxBinsEB; iNvtx++){ 
      sfEB[iEta][iNvtx] = dataEB[iEta][iNvtx]/mcEB[iEta][iNvtx];
      float sigmaDoDEB   = dataErrEB[iEta][iNvtx]/dataEB[iEta][iNvtx];
      float sigmaMCoMCEB = mcErrEB[iEta][iNvtx]/mcEB[iEta][iNvtx];
      sfErrTotEB[iEta][iNvtx] = sfEB[iEta][iNvtx]*sqrt( (sigmaDoDEB*sigmaDoDEB) + (sigmaMCoMCEB*sigmaMCoMCEB) );
      std::cout << sfEB[iEta][iNvtx] << " +/- " << sfErrTotEB[iEta][iNvtx] << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << "scale factors: EE" << std::endl;
  double sfEE[nEtaBins][nNvtxBinsEE];
  double sfErrTotEE[nEtaBins][nNvtxBinsEE];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iNvtx=0; iNvtx<nNvtxBinsEE; iNvtx++){ 
      sfEE[iEta][iNvtx] = dataEE[iEta][iNvtx]/mcEE[iEta][iNvtx];
      float sigmaDoDEE   = dataErrEE[iEta][iNvtx]/dataEE[iEta][iNvtx];
      float sigmaMCoMCEE = mcErrEE[iEta][iNvtx]/mcEE[iEta][iNvtx];
      sfErrTotEE[iEta][iNvtx] = sfEE[iEta][iNvtx]*sqrt( (sigmaDoDEE*sigmaDoDEE) + (sigmaMCoMCEE*sigmaMCoMCEE) );
      std::cout << sfEE[iEta][iNvtx] << " +/- " << sfErrTotEE[iEta][iNvtx] << std::endl;
    }
  }


  // Draw all canvases
  for(int ieta = 0; ieta<nEtaBins; ieta++){

    TString cname = "sfEff_";
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
    double *dataSliceEB    = new double[nNvtxBinsEB];
    double *dataSliceErrEB = new double[nNvtxBinsEB];
    double *mcSliceEB      = new double[nNvtxBinsEB];
    double *mcSliceErrEB   = new double[nNvtxBinsEB];
    double *sfSliceEB      = new double[nNvtxBinsEB];
    double *sfSliceErrEB   = new double[nNvtxBinsEB];
    for(int invtx = 0; invtx<nNvtxBinsEB; invtx++){
      dataSliceEB   [invtx] = dataEB     [ieta][invtx];
      dataSliceErrEB[invtx] = dataErrEB  [ieta][invtx];
      mcSliceEB     [invtx] = mcEB       [ieta][invtx];
      mcSliceErrEB  [invtx] = mcErrEB    [ieta][invtx];
      sfSliceEB     [invtx] = sfEB       [ieta][invtx];
      sfSliceErrEB  [invtx] = sfErrTotEB [ieta][invtx];
    }

    double *dataSliceEE    = new double[nNvtxBinsEE];
    double *dataSliceErrEE = new double[nNvtxBinsEE];
    double *mcSliceEE      = new double[nNvtxBinsEE];
    double *mcSliceErrEE   = new double[nNvtxBinsEE];
    double *sfSliceEE      = new double[nNvtxBinsEE];
    double *sfSliceErrEE   = new double[nNvtxBinsEE];
    for(int invtx = 0; invtx<nNvtxBinsEE; invtx++){
      dataSliceEE   [invtx] = dataEE     [ieta][invtx];
      dataSliceErrEE[invtx] = dataErrEE  [ieta][invtx];
      mcSliceEE     [invtx] = mcEE       [ieta][invtx];
      mcSliceErrEE  [invtx] = mcErrEE    [ieta][invtx];
      sfSliceEE     [invtx] = sfEE       [ieta][invtx];
      sfSliceErrEE  [invtx] = sfErrTotEE [ieta][invtx];
    }

    // Create and configure the graphs   
    TGraphErrors *grDataEB = new TGraphErrors(nNvtxBinsEB, nvtxBinCentersEB, dataSliceEB, nvtxBinHalfWidthEB, dataSliceErrEB);
    TGraphErrors *grMcEB   = new TGraphErrors(nNvtxBinsEB, nvtxBinCentersEB, mcSliceEB, nvtxBinHalfWidthEB, mcSliceErrEB);
    TGraphErrors *grSfEB   = new TGraphErrors(nNvtxBinsEB, nvtxBinCentersEB, sfSliceEB, nvtxBinHalfWidthEB, sfSliceErrEB);

    TGraphErrors *grDataEE = new TGraphErrors(nNvtxBinsEE, nvtxBinCentersEE, dataSliceEE, nvtxBinHalfWidthEE, dataSliceErrEE);
    TGraphErrors *grMcEE   = new TGraphErrors(nNvtxBinsEE, nvtxBinCentersEE, mcSliceEE, nvtxBinHalfWidthEE, mcSliceErrEE);
    TGraphErrors *grSfEE   = new TGraphErrors(nNvtxBinsEE, nvtxBinCentersEE, sfSliceEE, nvtxBinHalfWidthEE, sfSliceErrEE);
    
    grDataEB->SetLineColor(kBlack);
    grDataEB->SetMarkerColor(kBlack);
    grDataEB->SetMarkerStyle(20);
    grDataEB->SetMarkerSize(1.);
    grDataEE->SetLineColor(kBlack);
    grDataEE->SetMarkerColor(kBlack);
    grDataEE->SetMarkerStyle(20);
    grDataEE->SetMarkerSize(1.);

    int ci = TColor::GetColor("#99ccff");
    grMcEB->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grMcEB->SetLineColor(kGreen+4);
    grMcEB->SetMarkerStyle(22);
    grMcEB->SetMarkerColor(kGreen+4);
    grMcEB->SetMarkerSize(1.);

    ci = TColor::GetColor("#99ccff");
    grMcEE->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grMcEE->SetLineColor(kGreen+4);
    grMcEE->SetMarkerStyle(22);
    grMcEE->SetMarkerColor(kGreen+4);
    grMcEE->SetMarkerSize(1.);

    ci = TColor::GetColor("#99ccff");
    grSfEB->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grSfEB->SetLineColor(kGreen+4);
    grSfEB->SetMarkerStyle(20);
    grSfEB->SetMarkerColor(kGreen+4);
    grSfEB->SetMarkerSize(1.);
    //grSfEB->Fit("pol0","","",0,50);

    ci = TColor::GetColor("#99ccff");
    grSfEE->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grSfEE->SetLineColor(kGreen+4);
    grSfEE->SetMarkerStyle(20);
    grSfEE->SetMarkerColor(kGreen+4);
    grSfEE->SetMarkerSize(1.);
    //grSfEE->Fit("pol0","","",0,50);

    // Create and configure the dummy histograms on which to draw the graphs
    TH2F *h1 = new TH2F("dummy1","", 100, 0, 50, 100, 0.2, 1.1);
    h1->GetYaxis()->SetTitle("Efficiency");
    h1->SetStats(0);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetNdivisions(505);
    h1->GetXaxis()->SetDecimals();
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetYaxis()->SetTitleSize(0.05);
    TH2F *h2 = new TH2F("dummy2","", 100, 0, 50, 100, 0.8, 1.5);
    // h2->GetXaxis()->SetTitle("#vtx");
    h2->GetXaxis()->SetTitle("#rho");
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
    leg->AddEntry(grDataEB, "Data", "pl");
    leg->AddEntry(grMcEB, "Simulation DY", "pFlE");

    TLatex *latLumi = new TLatex(0, 1.15, lumiString);

    TLatex *latEtaEB = new TLatex(60.0, 0.5, etaLimitsStringArrayEB[ieta]);
    TLatex *latEtaEE = new TLatex(60.0, 0.5, etaLimitsStringArrayEE[ieta]);


    // --------------------------------------
    // EB
    // Draw the efficiencies
    pad1->cd();
    h1->Draw();
    grMcEB  ->Draw("2same");
    grMcEB  ->Draw("pZ,same");
    grDataEB->Draw("PEZ,same");
    leg->Draw("same");
    latEtaEB->Draw("same");
    latLumi->Draw("same");
    // Draw the scale factors
    pad2->cd();
    h2->Draw();
    grSfEB  ->Draw("2same");
    grSfEB  ->Draw("pEZ,same");
    // Save into a file
    TString fname = cname;
    fname += "_EB.pdf";
    c1->Print(fname);

    // --------------------------------------
    // EE
    // Draw the efficiencies
    pad1->cd();
    h1->Draw();
    grMcEE  ->Draw("2same");
    grMcEE  ->Draw("pZ,same");
    grDataEE->Draw("PEZ,same");
    leg->Draw("same");
    latEtaEE->Draw("same");
    latLumi->Draw("same");
    // Draw the scale factors
    pad2->cd();
    h2->Draw();
    grSfEE  ->Draw("2same");
    grSfEE  ->Draw("pEZ,same");
    // Save into a file
    fname = cname;
    fname += "_EE.pdf";
    c1->Print(fname);
  }

}



