#include "TString.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>

bool wantSyst = false;

const int nEtaBins = 1;

const TString lumiString = "CMS Preliminary, #sqrt{s}=13 TeV, #intLdt=0.801 fb^{-1}";

// EB
const int nPtBinsEB = 8;   
const double ptBinLimitsEB[nPtBinsEB+1]  = {20., 30., 40., 50., 60., 80., 110., 150., 200.,};
const double ptBinCentersEB[nPtBinsEB]   = {25., 35., 45., 55., 70., 95., 130., 175.};
const double ptBinHalfWidthEB[nPtBinsEB] = { 5.,  5.,  5.,  5., 10., 15.,  20.,  25.};
const TString etaLimitsStringArrayEB[nEtaBins] = { "0. < |#eta| < 1.4442" };

// EE
const int nPtBinsEE = 8;
const double ptBinLimitsEE[nPtBinsEE+1]  = {20., 30., 40., 50., 60., 80., 110., 150., 200.};
const double ptBinCentersEE[nPtBinsEE]   = {25., 35., 45., 55., 70., 95., 130., 175.};
const double ptBinHalfWidthEE[nPtBinsEE] = { 5.,  5.,  5.,  5., 10., 15.,  20.,  25.};
const TString etaLimitsStringArrayEE[nEtaBins] = { "1.566 < |#eta| < 2.5" };


// ----------------------------------------
// Data efficiencies and statistical errors
double dataEB[nEtaBins][nPtBinsEB] = {
  { 
    //7.90633e-01, 8.63725e-01, 8.88852e-01, 8.91585e-01, 8.94966e-01, 8.95945e-01, 8.85518e-01, 9.00705e-01   // ALL freeze
    //7.99912e-01, 8.69168e-01, 8.92712e-01, 8.95969e-01, 9.00560e-01, 9.02923e-01, 8.93823e-01, 8.99063e-01      // N-1 CHiso
    //9.01861e-01, 9.59721e-01, 9.76665e-01, 9.73847e-01, 9.73263e-01, 9.74641e-01, 9.68698e-01, 9.73443e-01      // N-1 PHiso
    //7.98860e-01, 8.67611e-01, 8.90312e-01, 8.92469e-01, 8.95844e-01, 8.96454e-01, 8.88051e-01, 9.02070e-01      // N-1 H/E
    8.11516e-01, 8.77539e-01, 9.00340e-01, 9.04088e-01, 9.07939e-01, 9.04841e-01, 8.98301e-01, 9.04669e-01      // N-1 sieie
  }
};

double dataEE[nEtaBins][nPtBinsEE] = {
  { 
    6.97114e-01, 7.58579e-01, 7.94582e-01, 7.85836e-01, 7.99429e-01, 8.57930e-01, 8.57685e-01, 8.86719e-01        // chiara: gli ultimi tre sono quelli vecchi
  }
};

// statistical only errors 
double dataErrStatEB[nEtaBins][nPtBinsEB] = {
  { 
    //4.00601e-03, 1.12189e-03, 7.68246e-04, 1.72324e-03, 2.96061e-03, 6.83292e-03, 1.26177e-02, 1.79148e-02
    //4.10345e-03, 1.11136e-03, 7.56770e-04, 1.69318e-03, 2.83079e-03, 6.76579e-03, 1.18051e-02, 2.10301e-02      // N-1 CHiso     
    //4.62149e-03, 9.85554e-04, 4.69735e-04, 8.65355e-04, 2.54654e-03, 5.73005e-03, 9.21217e-03, 1.54683e-02      // N-1 PHiso
    //3.95654e-03, 1.10462e-03, 7.65389e-04, 1.70874e-03, 2.95475e-03, 6.78395e-03, 1.26134e-02, 1.76768e-02      // N-1 H/E
    3.99652e-03, 1.08336e-03, 7.36002e-04, 1.63096e-03, 2.75232e-03, 6.41883e-03, 1.16012e-02, 1.73787e-02      // N-1 sieie    
  }
};

double dataErrStatEE[nEtaBins][nPtBinsEE] = {
  { 
    5.10971e-03, 2.17719e-03, 1.73005e-03, 4.42927e-03, 8.62274e-03, 1.89293e-02, 2.56017e-02, 3.64193e-02 
  }
};


// ----------------------------------------
// alternative fit changing the signal model and keeping nominal background
double dataSystSigEB[nEtaBins][nPtBinsEB] = {
  { 
    8.15170e-01, 8.60076e-01, 8.93184e-01, 8.96896e-01, 8.96205e-01, 9.11385e-01, 8.92919e-01, 9.04949e-01 //, 9.32198e-01, 9.30125e-01, 8.47152e-01 
  }
};

double dataSystSigEE[nEtaBins][nPtBinsEE] = {
  { 
    7.05637e-01, 7.57202e-01, 7.96266e-01, 8.06905e-01, 8.37237e-01, 8.60022e-01, 8.50684e-01, 8.91316e-01 //, 9.03238e-01 
  }
};


// ----------------------------------------
// alternative fit changing the background model and keeping nominal signal
double dataSystBackEB[nEtaBins][nPtBinsEB] = {
  { 
    8.23609e-01, 8.64064e-01, 8.90442e-01, 8.90713e-01, 8.94145e-01, 8.96122e-01, 9.05821e-01, 8.92899e-01 // , 9.22727e-01, 9.24031e-01, 8.68644e-01 
  }
};

double dataSystBackEE[nEtaBins][nPtBinsEE] = {
  { 
    6.91353e-01, 7.66270e-01, 7.98943e-01, 8.05622e-01, 8.21131e-01, 8.26827e-01, 8.97592e-01, 9.24601e-01 //, 9.17364e-01
  }
};
  

// ----------------------------------------
// MC efficiencies and errors - C&C
double mcEB[nEtaBins][nPtBinsEB] = {
  { 
    //0.846248, 0.888304, 0.911131, 0.913963, 0.914291, 0.912208, 0.91753,  0.916693   // FullSel
    //0.855629, 0.893929, 0.914799, 0.917692, 0.918364, 0.917077, 0.92329,  0.921731   // N-1 Chiso
    //0.941808, 0.972231, 0.983033, 0.982036, 0.980853, 0.978614, 0.978999, 0.978118   // N-1 Phiso
    //0.852861, 0.891288, 0.912293, 0.914747, 0.914996, 0.913228, 0.918136, 0.918424   // N-1 H/E
    0.864875, 0.899079, 0.919486, 0.923232, 0.923886, 0.921421, 0.925580, 0.924595   // N-1 sieie
  }
};

double mcEE[nEtaBins][nPtBinsEE] = {
  { 
    0.685745, 0.747946, 0.784424, 0.801722, 0.81622, 0.837517, 0.849242, 0.873572
  }
};

// statistical only errors - 218pb
double mcErrEB[nEtaBins][nPtBinsEB] = {
  { 
    0.00028343, 0.000149217, 0.000129587, 0.000274898, 0.000460165, 0.000945207, 0.00162603, 0.00266309
  }
};

double mcErrEE[nEtaBins][nPtBinsEE] = {
  { 
    0.00057193, 0.000345498, 0.000318749, 0.000676631, 0.00113528, 0.00227925, 0.00399321, 0.00638182
  }
};


void drawResults2016nm1(){

  // Syst error: take the max difference
  double dataSystErrEB[nEtaBins][nPtBinsEB];
  for (int ii=0; ii<nPtBinsEB; ii++ ) {
    if ( fabs(dataEB[0][ii]-dataSystSigEB[0][ii]) > fabs(dataEB[0][ii]-dataSystBackEB[0][ii]) ) dataSystErrEB[0][ii] = fabs(dataEB[0][ii]-dataSystSigEB[0][ii]);
    else dataSystErrEB[0][ii] = fabs(dataEB[0][ii]-dataSystBackEB[0][ii]);
  }

  double dataSystErrEE[nEtaBins][nPtBinsEE];
  for (int ii=0; ii<nPtBinsEE; ii++ ) {
    if ( fabs(dataEE[0][ii]-dataSystSigEE[0][ii]) > fabs(dataEE[0][ii]-dataSystBackEE[0][ii]) ) dataSystErrEE[0][ii] = fabs(dataEE[0][ii]-dataSystSigEE[0][ii]);
    else dataSystErrEE[0][ii] = fabs(dataEE[0][ii]-dataSystBackEE[0][ii]);
  }

  // Tot error: stat + syst
  double dataErrEB[nEtaBins][nPtBinsEB];
  double dataErrEE[nEtaBins][nPtBinsEE];

  if (wantSyst) {
    std::cout << "systematics added" << std::endl;
    for (int ii=0; ii<nPtBinsEB; ii++ ) {
      dataErrEB[0][ii] = sqrt( dataSystErrEB[0][ii]*dataSystErrEB[0][ii] + dataErrStatEB[0][ii]*dataErrStatEB[0][ii] );
    }
    for (int ii=0; ii<nPtBinsEE; ii++ ) {
      dataErrEE[0][ii] = sqrt( dataSystErrEE[0][ii]*dataSystErrEE[0][ii] + dataErrStatEE[0][ii]*dataErrStatEE[0][ii] );
    }
  } else {
    std::cout << "only statistical error" << std::endl;
    for (int ii=0; ii<nPtBinsEB; ii++ ) {
      dataErrEB[0][ii] = dataErrStatEB[0][ii];
    }
    for (int ii=0; ii<nPtBinsEE; ii++ ) {
      dataErrEE[0][ii] = dataErrStatEE[0][ii];
    }
  }

  std::cout << "================================" << std::endl;
  std::cout << "EB" << std::endl;
  for (int ii=0; ii<nPtBinsEB; ii++ ) 
    std::cout << ii << ", nominal = " << dataEB[0][ii] << ", forSigSyst = " << dataSystSigEB[0][ii] << ", forBkgSyst = " << dataSystBackEB[0][ii] 
	 << ", statErr = " << dataErrStatEB[0][ii] << ", systErr = " <<dataSystErrEB[0][ii] << std::endl; 
  std::cout << "================================" << std::endl;
  std::cout << "EE" << std::endl;
  for (int ii=0; ii<nPtBinsEE; ii++ ) 
    std::cout << ii << ", nominal = " << dataEE[0][ii] << ", forSigSyst = " << dataSystSigEE[0][ii] << ", forBkgSyst = " << dataSystBackEE[0][ii] 
	 << ", statErr = " << dataErrStatEE[0][ii] << ", systErr = " <<dataSystErrEE[0][ii] << std::endl; 
  std::cout << "================================" << std::endl;
  

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << "scale factors: EB" << std::endl;
  // Scale factors and errors
  double sfEB[nEtaBins][nPtBinsEB];
  double sfErrTotEB[nEtaBins][nPtBinsEB];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iPt=0; iPt<nPtBinsEB; iPt++){ 
      sfEB[iEta][iPt] = dataEB[iEta][iPt]/mcEB[iEta][iPt];
      float sigmaDoDEB   = dataErrEB[iEta][iPt]/dataEB[iEta][iPt];
      float sigmaMCoMCEB = mcErrEB[iEta][iPt]/mcEB[iEta][iPt];
      sfErrTotEB[iEta][iPt] = sfEB[iEta][iPt]*sqrt( (sigmaDoDEB*sigmaDoDEB) + (sigmaMCoMCEB*sigmaMCoMCEB) );
      std::cout << sfEB[iEta][iPt] << " +/- " << sfErrTotEB[iEta][iPt] << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << "scale factors: EE" << std::endl;
  double sfEE[nEtaBins][nPtBinsEE];
  double sfErrTotEE[nEtaBins][nPtBinsEE];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iPt=0; iPt<nPtBinsEE; iPt++){ 
      sfEE[iEta][iPt] = dataEE[iEta][iPt]/mcEE[iEta][iPt];
      float sigmaDoDEE   = dataErrEE[iEta][iPt]/dataEE[iEta][iPt];
      float sigmaMCoMCEE = mcErrEE[iEta][iPt]/mcEE[iEta][iPt];
      sfErrTotEE[iEta][iPt] = sfEE[iEta][iPt]*sqrt( (sigmaDoDEE*sigmaDoDEE) + (sigmaMCoMCEE*sigmaMCoMCEE) );
      std::cout << sfEE[iEta][iPt] << " +/- " << sfErrTotEE[iEta][iPt] << std::endl;
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
    double *dataSliceEB    = new double[nPtBinsEB];
    double *dataSliceErrEB = new double[nPtBinsEB];
    double *mcSliceEB      = new double[nPtBinsEB];
    double *mcSliceErrEB   = new double[nPtBinsEB];
    double *sfSliceEB      = new double[nPtBinsEB];
    double *sfSliceErrEB   = new double[nPtBinsEB];
    for(int ipt = 0; ipt<nPtBinsEB; ipt++){
      dataSliceEB   [ipt] = dataEB     [ieta][ipt];
      dataSliceErrEB[ipt] = dataErrEB  [ieta][ipt];
      mcSliceEB     [ipt] = mcEB       [ieta][ipt];
      mcSliceErrEB  [ipt] = mcErrEB    [ieta][ipt];
      sfSliceEB     [ipt] = sfEB       [ieta][ipt];
      sfSliceErrEB  [ipt] = sfErrTotEB [ieta][ipt];
    }

    double *dataSliceEE    = new double[nPtBinsEE];
    double *dataSliceErrEE = new double[nPtBinsEE];
    double *mcSliceEE      = new double[nPtBinsEE];
    double *mcSliceErrEE   = new double[nPtBinsEE];
    double *sfSliceEE      = new double[nPtBinsEE];
    double *sfSliceErrEE   = new double[nPtBinsEE];
    for(int ipt = 0; ipt<nPtBinsEE; ipt++){
      dataSliceEE   [ipt] = dataEE     [ieta][ipt];
      dataSliceErrEE[ipt] = dataErrEE  [ieta][ipt];
      mcSliceEE     [ipt] = mcEE       [ieta][ipt];
      mcSliceErrEE  [ipt] = mcErrEE    [ieta][ipt];
      sfSliceEE     [ipt] = sfEE       [ieta][ipt];
      sfSliceErrEE  [ipt] = sfErrTotEE [ieta][ipt];
    }

    // Create and configure the graphs   
    TGraphErrors *grDataEB = new TGraphErrors(nPtBinsEB, ptBinCentersEB, dataSliceEB, ptBinHalfWidthEB, dataSliceErrEB);
    TGraphErrors *grMcEB   = new TGraphErrors(nPtBinsEB, ptBinCentersEB, mcSliceEB, ptBinHalfWidthEB, mcSliceErrEB);
    TGraphErrors *grSfEB   = new TGraphErrors(nPtBinsEB, ptBinCentersEB, sfSliceEB, ptBinHalfWidthEB, sfSliceErrEB);

    TGraphErrors *grDataEE = new TGraphErrors(nPtBinsEE, ptBinCentersEE, dataSliceEE, ptBinHalfWidthEE, dataSliceErrEE);
    TGraphErrors *grMcEE   = new TGraphErrors(nPtBinsEE, ptBinCentersEE, mcSliceEE, ptBinHalfWidthEE, mcSliceErrEE);
    TGraphErrors *grSfEE   = new TGraphErrors(nPtBinsEE, ptBinCentersEE, sfSliceEE, ptBinHalfWidthEE, sfSliceErrEE);
    
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

    ci = TColor::GetColor("#99ccff");
    grSfEE->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grSfEE->SetLineColor(kGreen+4);
    grSfEE->SetMarkerStyle(20);
    grSfEE->SetMarkerColor(kGreen+4);
    grSfEE->SetMarkerSize(1.);

    // Create and configure the dummy histograms on which to draw the graphs
    // TH2F *h1 = new TH2F("dummy1","", 100, 0, 500, 100, 0.6, 1.1);
    TH2F *h1 = new TH2F("dummy1","", 100, 0, 200, 100, 0.6, 1.1);
    h1->GetYaxis()->SetTitle("Efficiency");
    h1->SetStats(0);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetNdivisions(505);
    h1->GetXaxis()->SetDecimals();
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetYaxis()->SetTitleSize(0.05);
    //TH2F *h2 = new TH2F("dummy2","", 100, 0, 500, 100, 0.8, 1.2);
    TH2F *h2 = new TH2F("dummy2","", 100, 0, 200, 100, 0.8, 1.2);
    h2->GetXaxis()->SetTitle("p_{T} [GeV]");
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



