#include "TString.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TVirtualFitter.h"
#include <iostream>

bool wantSyst = true;

const int nEtaBins = 1;

const TString lumiString = "CMS Preliminary, #sqrt{s}=13 TeV, #intLdt=4.0 fb^{-1}";

// EB
const int nPtBinsEB = 10;   
const double ptBinLimitsEB[nPtBinsEB+1]  = {20., 30., 40., 50., 60., 80., 110., 150., 200.,270.,350.};
const double ptBinCentersEB[nPtBinsEB]   = {25., 35., 45., 55., 70., 95., 130., 175., 235., 310.};
const double ptBinHalfWidthEB[nPtBinsEB] = { 5.,  5.,  5.,  5., 10., 15.,  20.,  25., 35.,  40.};
const TString etaLimitsStringArrayEB[nEtaBins] = { "0. < |#eta| < 1.4442" };
//const int nPtBinsEB = 9;   
//const double ptBinLimitsEB[nPtBinsEB+1]  = {20., 30., 40., 50., 60., 80., 110., 150., 200.,270.};
//const double ptBinCentersEB[nPtBinsEB]   = {25., 35., 45., 55., 70., 95., 130., 175., 235.};
//const double ptBinHalfWidthEB[nPtBinsEB] = { 5.,  5.,  5.,  5., 10., 15.,  20.,  25., 35.};
//const TString etaLimitsStringArrayEB[nEtaBins] = { "0. < |#eta| < 1.4442" };

// EE
//const int nPtBinsEE = 8;
//const double ptBinLimitsEE[nPtBinsEE+1]  = {20., 30., 40., 50., 60., 80., 110., 150., 200.};
//const double ptBinCentersEE[nPtBinsEE]   = {25., 35., 45., 55., 70., 95., 130., 175.};
//const double ptBinHalfWidthEE[nPtBinsEE] = { 5.,  5.,  5.,  5., 10., 15.,  20.,  25.};
//const TString etaLimitsStringArrayEE[nEtaBins] = { "1.566 < |#eta| < 2.5" };
const int nPtBinsEE = 9;
const double ptBinLimitsEE[nPtBinsEE+1]  = {20., 30., 40., 50., 60., 80., 110., 150., 200.,350.};
const double ptBinCentersEE[nPtBinsEE]   = {25., 35., 45., 55., 70., 95., 130., 175., 275.};
const double ptBinHalfWidthEE[nPtBinsEE] = { 5.,  5.,  5.,  5., 10., 15.,  20.,  25., 75.};
const TString etaLimitsStringArrayEE[nEtaBins] = { "1.566 < |#eta| < 2.5" };


// ----------------------------------------
// Data efficiencies and statistical errors
double dataEB[nEtaBins][nPtBinsEB] = {
  { 
    7.91251e-01, 8.52342e-01, 8.82546e-01, 8.83756e-01, 8.79436e-01, 8.83524e-01, 8.85466e-01, 8.98044e-01, 8.86877e-01, 8.52180e-01
  }
};

double dataEE[nEtaBins][nPtBinsEE] = {
  { 
    // |eta|<2.5
    // 6.93590e-01, 7.66276e-01, 8.02283e-01, 8.15095e-01, 8.28115e-01, 8.41396e-01, 8.30929e-01, 8.44886e-01, 8.96420e-01 

    // |eta|<2.1
    6.67763e-01, 7.47609e-01, 7.91513e-01, 8.09517e-01, 8.18605e-01, 8.39861e-01, 8.35895e-01, 8.36613e-01, 8.98649e-01
  }
};

// statistical only errors 
double dataErrStatEB[nEtaBins][nPtBinsEB] = {
  { 
    1.77659e-03, 7.32544e-04, 4.30639e-04, 1.03912e-03, 1.92876e-03, 3.91642e-03, 6.39313e-03, 1.00293e-02, 1.45479e-02, 2.48560e-02
  }
};

double dataErrStatEE[nEtaBins][nPtBinsEE] = {
  { 
    // |eta|<2.5    
    // 2.26404e-03, 1.14846e-03, 1.27566e-04, 2.54440e-03, 4.39143e-03, 8.51082e-03, 1.30638e-02, 1.87204e-02, 2.32021e-02

    // |eta|<2.1  
    2.89224e-03, 1.42109e-03, 1.27520e-03, 2.79981e-03, 5.34937e-03, 9.76875e-03, 1.49404e-02, 2.40042e-02, 2.45539e-02
  }
};


// ----------------------------------------
// alternative fit changing the signal model and keeping nominal background
double dataSystSigEB[nEtaBins][nPtBinsEB] = {
  { 
    8.20592e-01, 8.50285e-01, 8.81259e-01, 8.86457e-01, 8.89836e-01, 8.98114e-01, 8.89983e-01, 9.02391e-01, 8.75210e-01, 8.25859e-01
  }
};

double dataSystSigEE[nEtaBins][nPtBinsEE] = {
  { 
    // |eta|<2.1  
    6.83923e-01, 7.54630e-01, 7.88666e-01, 8.10450e-01, 8.26529e-01, 8.49273e-01, 8.41623e-01, 8.18942e-01, 9.00141e-01
  }
};


// ----------------------------------------
// alternative fit changing the background model and keeping nominal signal
double dataSystBackEB[nEtaBins][nPtBinsEB] = {
  { 
    8.01937e-01, 8.52328e-01, 8.82587e-01, 8.84032e-01, 8.79999e-01, 8.83986e-01, 8.85510e-01, 9.04938e-01, 8.91095e-01, 8.52497e-01
  }
};

double dataSystBackEE[nEtaBins][nPtBinsEE] = {
  { 
    // |eta|<2.1   
    6.68713e-01, 7.47641e-01, 7.91769e-01, 8.09870e-01, 8.20948e-01, 8.46312e-01, 8.35896e-01, 8.37405e-01, 9.03088e-01
  }
};


// ----------------------------------------
// MC efficiencies and errors - C&C
double mcEB[nEtaBins][nPtBinsEB] = {
  { 
    0.832083, 0.880983, 0.908311, 0.910573, 0.910272, 0.910244, 0.914513, 0.915982, 0.928903, 0.926116
  }
};

double mcEE[nEtaBins][nPtBinsEE] = {
  { 
    // |eta|<2.5   
    // 0.68196, 0.74986, 0.788325, 0.805504, 0.818779, 0.838306, 0.850299, 0.882256, 0.894807

    // |eta|<2.1
    0.68177,  0.7528,   0.795943, 0.813491, 0.828562, 0.851386, 0.864961, 0.891915, 0.914476
  }
};

// statistical only errors - 218pb
double mcErrEB[nEtaBins][nPtBinsEB] = {
  { 
    0.00028343, 0.000149217, 0.000129587, 0.000274898, 0.000460165, 0.000945207, 0.00162603, 0.00266309, 0.0041, 0.0074   // ultimi due ottenuti dividendo moriond per 1.7, che e' il rapporto tra i due nei bin precedenti
  }
};

double mcErrEE[nEtaBins][nPtBinsEE] = {
  { 
    0.00057193, 0.000345498, 0.000318749, 0.000676631, 0.00113528, 0.00227925, 0.00399321, 0.00638182, 0.008    // ultimi ottenuti dividendo moriond per 1.7. Uso gli stessi anche per |eta|<2.1
  }
};


void drawResults(){

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

    // Create an d configure the graphs   
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
    grSfEB->Fit("pol0","","",20,350);

    ci = TColor::GetColor("#99ccff");
    grSfEE->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grSfEE->SetLineColor(kGreen+4);
    grSfEE->SetMarkerStyle(20);
    grSfEE->SetMarkerColor(kGreen+4);
    grSfEE->SetMarkerSize(1.);
    grSfEE->Fit("pol0","","",20,350);

    // Create and configure the dummy histograms on which to draw the graphs
    // TH2F *h1 = new TH2F("dummy1","", 100, 0, 500, 100, 0.6, 1.1);
    TH2F *h1 = new TH2F("dummy1","", 100, 0, 350, 100, 0.6, 1.1);
    h1->GetYaxis()->SetTitle("Efficiency");
    h1->SetStats(0);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetNdivisions(505);
    h1->GetXaxis()->SetDecimals();
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetYaxis()->SetTitleSize(0.05);
    //TH2F *h2 = new TH2F("dummy2","", 100, 0, 500, 100, 0.8, 1.2);
    TH2F *h2 = new TH2F("dummy2","", 100, 0, 350, 100, 0.8, 1.2);
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


  // Now efficiency error study - from Tommaso

  // barrel SFs ----------------------------------------------
  TGraphErrors *gEB = new TGraphErrors();
  for (int iPt=0; iPt<nPtBinsEB; iPt++){   
    gEB->SetPoint(iPt,ptBinCentersEB[iPt],sfEB[0][iPt]);
    gEB->SetPointError(iPt,ptBinHalfWidthEB[iPt],sfErrTotEB[0][iPt]);
  }

  TCanvas *cSFEB = new TCanvas("cSFEB","cSFEB");
  cSFEB->cd();
  gPad->SetGrid();
  gEB->SetMarkerStyle(20);
  int ci = TColor::GetColor("#AAFF55");
  gEB->SetFillColor(ci);
  TH1F *hEB = gPad->DrawFrame(0.,0.8,360.,1.2); 
  hEB->GetXaxis()->SetTitle("p_{T} [GeV]");
  hEB->GetYaxis()->SetTitle("Scale Factor"); 
  hEB->Draw();
  gEB->Draw("P");

  // pol0 fit
  TF1 *userEB = new TF1("userEB","[0]+[1]*x",40.,1000.);
  userEB->SetParameters(1.,0); 
  userEB->FixParameter(1.,0); 
  gEB->Fit("userEB","R");
  TH1F *hIntEB = new TH1F("hIntEB","",100,40.,1000);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hIntEB,0.68);
  hIntEB -> SetStats(kFALSE);
  hIntEB -> SetFillStyle(3004);
  hIntEB -> SetFillColor(kBlue+1);
  hIntEB ->Draw("E3same");

  // pol1 fit
  TF1 *user1EB = new TF1("user1EB","[0]+[1]*x",40.,1000.);
  user1EB->SetParameters(1.,0.); 
  gEB->Fit("user1EB","R");
  TH1F *hInt1EB = new TH1F("hInt1EB","",100,40.,1000);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hInt1EB,0.68);
  hInt1EB -> SetStats(kFALSE);
  hInt1EB -> SetFillStyle(3005);
  hInt1EB -> SetFillColor(kRed+1);
  hInt1EB ->Draw("E3same");

  TLegend *tlEB = new TLegend(0.17,0.18,0.65,0.3);
  tlEB->SetFillColor(0);
  tlEB->SetTextFont(42);
  tlEB->AddEntry(hIntEB, "Pol0 - plateau region: 68% CL","F");
  tlEB->AddEntry(hInt1EB,"Pol1 - plateau region: 68% CL","F");
  tlEB->Draw();

  cSFEB->SaveAs("SFuncertaintyEB.png");


  // endcpas SFs ----------------------------------------------
  TGraphErrors *gEE = new TGraphErrors();
  for (int iPt=0; iPt<nPtBinsEE; iPt++){   
    gEE->SetPoint(iPt,ptBinCentersEE[iPt],sfEE[0][iPt]);
    gEE->SetPointError(iPt,ptBinHalfWidthEE[iPt],sfErrTotEE[0][iPt]);
  }

  TCanvas *cSFEE = new TCanvas("cSFEE","cSFEE");
  cSFEE->cd();
  gPad->SetGrid();
  gEE->SetMarkerStyle(20);
  ci = TColor::GetColor("#AAFF55");
  gEE->SetFillColor(ci);
  TH1F *hEE = gPad->DrawFrame(0.,0.8,360.,1.2); 
  hEE->GetXaxis()->SetTitle("p_{T} [GeV]");
  hEE->GetYaxis()->SetTitle("Scale Factor"); 
  hEE->Draw();
  gEE->Draw("P");

  // pol0 fit
  TF1 *userEE = new TF1("userEE","[0]+[1]*x",40.,1000.);
  userEE->SetParameters(1.,0); 
  userEE->FixParameter(1.,0); 
  // gEE->Fit("userEE","R");
  gEE->Fit("userEE","","",75.,1000.);
  TH1F *hIntEE = new TH1F("hIntEE","",100,40.,1000);            
  // TH1F *hIntEE = new TH1F("hIntEE","",100,75.,1000);                    // chiara: solo per il plot dello SF    
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hIntEE,0.68);
  hIntEE -> SetStats(kFALSE);
  hIntEE -> SetFillStyle(3004);
  hIntEE -> SetFillColor(kBlue+1);
  hIntEE ->Draw("E3same");

  // pol1 fit
  TF1 *user1EE = new TF1("user1EE","[0]+[1]*x",40.,1000.);  
  user1EE->SetParameters(1.,0.); 
  //gEE->Fit("user1EE","R");
  gEE->Fit("user1EE","","",75.,1000.);
  TH1F *hInt1EE = new TH1F("hInt1EE","",100,40.,1000);          
  // TH1F *hInt1EE = new TH1F("hInt1EE","",100,75.,1000);                   // chiara: solo per il plot dello SF  
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hInt1EE,0.68);
  hInt1EE -> SetStats(kFALSE);
  hInt1EE -> SetFillStyle(3005);
  hInt1EE -> SetFillColor(kRed+1);
  hInt1EE ->Draw("E3same");

  TLegend *tlEE = new TLegend(0.17,0.18,0.65,0.3);
  tlEE->SetFillColor(0);
  tlEE->SetTextFont(42);
  tlEE->AddEntry(hIntEE, "Pol0 - plateau region: 68% CL","F");
  tlEE->AddEntry(hInt1EE,"Pol1 - plateau region: 68% CL","F");
  tlEE->Draw();

  cSFEE->SaveAs("SFuncertaintyEE.png");



  // EBEB category : 
  TCanvas *c1EBEB = new TCanvas("c1EBEB","c1EBEB");
  c1EBEB->cd(); 
  gPad->SetGridy(); 
  gPad->SetGridx(); 
  TH1F *hIntEEBEB = new TH1F("hIntEEBEB","",100,80.,2000);
  TH1F *hIntTEBEB = new TH1F("hIntTEBEB","",100,80.,2000);
  for (int i=1; i<hInt1EB->GetNbinsX(); i++) {
    hIntEEBEB->SetBinContent(i,0); 
    float xx = hInt1EB -> GetBinCenter(i); 
    float tt = 2. * sqrt(hInt1EB->GetBinError(i)*hInt1EB->GetBinError(i) + 
			 pow((user1EB->Eval(xx)-userEB->Eval(xx)),2));
    float ee = 2. * hInt1EB->GetBinError(i); 
    hIntEEBEB->SetBinError(i,100*ee); 
    hIntTEBEB->SetBinError(i,100*tt); 

    if (hIntEEBEB->GetBinLowEdge(i)<750 && (hIntEEBEB->GetBinLowEdge(i)+hIntEEBEB->GetBinWidth(i))>750) 
      cout << "this is the 750 bin: in EBEB, " << hIntEEBEB->GetBinError(i) << " " << hIntTEBEB->GetBinError(i) << endl; 
  }

  TH1F *hMasEBEB = (TH1F*)gPad->DrawFrame(500.,-35.,1900.,35.);
  hMasEBEB->SetTitle("Category: EBEB"); 
  hMasEBEB -> GetXaxis() -> SetTitle("m#lower[-0.5]{_{#gamma#gamma}} (GeV)"); 
  hMasEBEB -> GetYaxis() -> SetTitle("Efficiency uncertainty (%)"); 
  hIntTEBEB->SetFillColor(kBlue+1);
  hIntTEBEB->SetFillStyle(3001);
  hIntTEBEB->Draw("E3same"); 
  hIntEEBEB->SetFillColor(kRed+1);
  hIntEEBEB->SetFillStyle(3005);
  hIntEEBEB->Draw("E3same"); 

  TLegend *tsEBEB = new TLegend(0.17,0.18,0.83,0.3);
  tsEBEB->SetFillColor(0);
  tsEBEB->SetTextFont(42);
  tsEBEB->AddEntry(hIntEEBEB,"EBEB: Pol1 68% CL","F");
  tsEBEB->AddEntry(hIntTEBEB,"EBEB: Pol1 68% CL #oplus #Delta(Pol1-Pol0)","F");
  tsEBEB->Draw();
  
  c1EBEB->SaveAs("MassUncertaintyEBEB.png");


  // EBEE category : 
  TCanvas *c1EBEE = new TCanvas("c1EBEE","c1EBEE");
  c1EBEE->cd(); 
  gPad->SetGridy(); 
  gPad->SetGridx(); 
  TH1F *hIntEEBEE = new TH1F("hIntEEBEE","",100,80.,2000);
  TH1F *hIntTEBEE = new TH1F("hIntTEBEE","",100,80.,2000);
  for (int i=1; i<hInt1EB->GetNbinsX(); i++) {
    hIntEEBEE->SetBinContent(i,0); 
    float xx = hInt1EB -> GetBinCenter(i); 
    float ttEB = sqrt(hInt1EB->GetBinError(i)*hInt1EB->GetBinError(i) + pow((user1EB->Eval(xx)-userEB->Eval(xx)),2));
    float ttEE = sqrt(hInt1EE->GetBinError(i)*hInt1EE->GetBinError(i) + pow((user1EE->Eval(xx)-userEE->Eval(xx)),2));
    float tt = sqrt( ttEB*ttEB + ttEE*ttEE);  
    float eeEB = hInt1EB->GetBinError(i); 
    float eeEE = hInt1EB->GetBinError(i); 
    float ee   = sqrt( eeEB*eeEB + eeEE*eeEE );
    hIntEEBEE->SetBinError(i,100*ee); 
    hIntTEBEE->SetBinError(i,100*tt); 

    if (hIntEEBEE->GetBinLowEdge(i)<750 && (hIntEEBEE->GetBinLowEdge(i)+hIntEEBEE->GetBinWidth(i))>750) 
      cout << "this is the 750 bin: in EBEE, " << hIntEEBEE->GetBinError(i) << " " << hIntTEBEE->GetBinError(i) << endl; 
  }

  TH1F *hMasEBEE = (TH1F*)gPad->DrawFrame(500.,-35.,1900.,35.);
  hMasEBEE->SetTitle("Category: EBEE"); 
  hMasEBEE -> GetXaxis() -> SetTitle("m#lower[-0.5]{_{#gamma#gamma}} (GeV)"); 
  hMasEBEE -> GetYaxis() -> SetTitle("Efficiency uncertainty (%)"); 
  hIntTEBEE->SetFillColor(kBlue+1);
  hIntTEBEE->SetFillStyle(3001);
  hIntTEBEE->Draw("E3same"); 
  hIntEEBEE->SetFillColor(kRed+1);
  hIntEEBEE->SetFillStyle(3005);
  hIntEEBEE->Draw("E3same"); 

  TLegend *tsEBEE = new TLegend(0.17,0.18,0.83,0.3);
  tsEBEE->SetFillColor(0);
  tsEBEE->SetTextFont(42);
  tsEBEE->AddEntry(hIntEEBEE,"EBEE: Pol1 68% CL","F");
  tsEBEE->AddEntry(hIntTEBEE,"EBEE: Pol1 68% CL #oplus #Delta(Pol1-Pol0)","F");
  tsEBEE->Draw();
  
  c1EBEE->SaveAs("MassUncertaintyEBEE.png");
}



