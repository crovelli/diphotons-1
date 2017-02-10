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
const int nPtBinsEB = 10;   
const double ptBinLimitsEB[nPtBinsEB+1]  = {20., 30., 40., 50., 60., 80., 110., 150., 200.,270.,350.};
const double ptBinCentersEB[nPtBinsEB]   = {25., 35., 45., 55., 70., 95., 130., 175., 235., 310.};
const double ptBinHalfWidthEB[nPtBinsEB] = { 5.,  5.,  5.,  5., 10., 15.,  20.,  25., 35.,  40.};
const TString etaLimitsStringArrayEB[nEtaBins] = { "0. < |#eta| < 1.4442" };

// ----------------------------------------
// Data efficiencies and statistical errors
double dataEB[nEtaBins][nPtBinsEB] = {
  { 
    // nominal
    // 7.91251e-01, 8.52342e-01, 8.82546e-01, 8.83756e-01, 8.79436e-01, 8.83524e-01, 8.85466e-01, 8.98044e-01, 8.86877e-01, 8.52180e-01

    // cut=3
    // 8.13670e-01, 8.75048e-01, 9.04549e-01, 9.02997e-01, 8.98248e-01, 9.00984e-01, 8.98954e-01, 9.08336e-01, 8.98284e-01, 8.68432e-01

    // cut=3.5
    // 8.41916e-01, 9.03530e-01, 9.31272e-01, 9.28072e-01, 9.22797e-01, 9.25544e-01, 9.17558e-01, 9.25379e-01, 9.16782e-01, 8.75255e-01

    // N-1
    9.09851e-01, 9.55844e-01, 9.71967e-01, 9.67031e-01, 9.63853e-01, 9.70270e-01, 9.67097e-01, 9.69841e-01, 9.50438e-01, 9.49781e-01
  }
};

// statistical only errors 
double dataErrStatEB[nEtaBins][nPtBinsEB] = {
  { 
    // nominal
    // 1.77659e-03, 7.32544e-04, 4.30639e-04, 1.03912e-03, 1.92876e-03, 3.91642e-03, 6.39313e-03, 1.00293e-02, 1.45479e-02, 2.48560e-02

    // cut=3  
    // 1.83205e-03, 6.00164e-04, 3.67356e-04, 8.20071e-04, 1.47789e-03, 3.84336e-03, 6.16291e-03, 1.00e-02, 1.40266e-02, 2.37675e-02  

    // cut=3.5 
    // 1.93717e-03, 5.91940e-04, 3.20663e-04, 8.98274e-04, 1.76837e-03, 3.55081e-03, 5.91162e-03, 7.73590e-03, 1.32638e-02, 2.33856e-02

    // N-1
    2.14891e-03, 4.46852e-04, 3.40791e-04, 7.10374e-04, 1.81910e-03, 3.25267e-03, 5.28233e-03, 2.97800e-03,  1.09814e-02, 1.79140e-02
  }
};

// ----------------------------------------
// MC efficiencies and errors - C&C
double mcEB[nEtaBins][nPtBinsEB] = {
  { 
    // nominal
    // 0.832083, 0.880983, 0.908311, 0.910573, 0.910272, 0.910244, 0.914513, 0.915982, 0.928903, 0.926116

    // cut=3
    // 0.851694, 0.899912, 0.925754, 0.926702, 0.925317, 0.923841, 0.926248, 0.926582, 0.936605, 0.934106

    // cut=3.5     
    // 0.876195, 0.923383, 0.947264, 0.946783, 0.944082, 0.941452, 0.941045, 0.938949, 0.948843, 0.945177

    // N-1
    0.930832, 0.966418, 0.979758, 0.9782, 0.976215, 0.975835, 0.976102, 0.975371, 0.979472, 0.983065
  }
};

// statistical only errors - 218pb
double mcErrEB[nEtaBins][nPtBinsEB] = {
  { 
    // Uso nominal ovunque. Ultimi due ottenuti dividendo moriond per 1.7, che e' il rapporto tra i due nei bin precedenti
    0.00028343, 0.000149217, 0.000129587, 0.000274898, 0.000460165, 0.000945207, 0.00162603, 0.00266309, 0.0041, 0.0074   
  }
};

void drawResults_phIsoTuning(){

  // Tot error: stat + syst
  double dataErrEB[nEtaBins][nPtBinsEB];

  std::cout << "only statistical error" << std::endl;
  for (int ii=0; ii<nPtBinsEB; ii++ ) {
    dataErrEB[0][ii] = dataErrStatEB[0][ii];
  }

  std::cout << "================================" << std::endl;
  std::cout << "EB" << std::endl;
  for (int ii=0; ii<nPtBinsEB; ii++ ) 
    std::cout << ii << ", nominal = " << dataEB[0][ii] << ", statErr = " << dataErrStatEB[0][ii] << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << "scale factors: EB" << std::endl;
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

    // Create and configure the graphs   
    TGraphErrors *grDataEB = new TGraphErrors(nPtBinsEB, ptBinCentersEB, dataSliceEB, ptBinHalfWidthEB, dataSliceErrEB);
    TGraphErrors *grMcEB   = new TGraphErrors(nPtBinsEB, ptBinCentersEB, mcSliceEB, ptBinHalfWidthEB, mcSliceErrEB);
    TGraphErrors *grSfEB   = new TGraphErrors(nPtBinsEB, ptBinCentersEB, sfSliceEB, ptBinHalfWidthEB, sfSliceErrEB);
    
    grDataEB->SetLineColor(kBlack);
    grDataEB->SetMarkerColor(kBlack);
    grDataEB->SetMarkerStyle(20);
    grDataEB->SetMarkerSize(1.);

    int ci = TColor::GetColor("#99ccff");
    grMcEB->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grMcEB->SetLineColor(kGreen+4);
    grMcEB->SetMarkerStyle(22);
    grMcEB->SetMarkerColor(kGreen+4);
    grMcEB->SetMarkerSize(1.);

    ci = TColor::GetColor("#99ccff");
    grSfEB->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grSfEB->SetLineColor(kGreen+4);
    grSfEB->SetMarkerStyle(20);
    grSfEB->SetMarkerColor(kGreen+4);
    grSfEB->SetMarkerSize(1.);
    grSfEB->Fit("pol0","","",20,350);

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
  }

}



