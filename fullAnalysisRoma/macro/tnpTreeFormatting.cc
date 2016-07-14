#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <TRandom.h>

using namespace std;

float scaleCorr(int run, float r9, float sceta);
float smearings(float r9, float sceta);
bool isPhoIsoOk(float sceta, float r9, float corrphoiso);
float gammaPtOnlyCorrPhIso(float sceta, float pt, float phoiso);
int effectiveAreaGammaRegion(float sceta);

void tnpTreeFormat(const char* filename, float lumiForW) {

  // initialize the random number
  TRandom myRandom(12345);

  cout << "Formatting " << filename << endl;  

  int increm=0;

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;
  TH1F  *h_sumW = 0;

  fileOrig = TFile::Open(TString("/cmsrm/pc28_2/crovelli/data/Exo/topup7/")+TString(filename));
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("tnpAna/TaPtree");
    h_sumW = (TH1F*)fileOrig->Get("tnpAna/h_sumW");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  fileOrig->cd();
  if (!treeOrig) {
    cout << "Tree tnpAna/TaPTree not existing !" << endl; 
    return;    
  }

  treeOrig->SetMakeClass(0);
  cout << "TreeOrig->Size = "<< treeOrig->GetEntries() << endl;

  // number of entries saved in the first tree
  int nentriesOrig = treeOrig->GetEntries();   

  // Tree for the final format
  TFile *fileNew = TFile::Open(TString("/tmp/crovelli/Formatted_")+TString(filename),"recreate");
  fileNew->ls();
  fileNew->cd();
  TDirectory *myDir = (TDirectory*)fileNew->mkdir("tnpAna");
  myDir->cd();
  TTree *treeNew = new TTree("TaPTree","reduced tree for T&P");
  treeNew->SetAutoSave(-99999999999);
  treeNew->SetAutoFlush(-99999999999);
  
  std::vector<TTree*> trees; 
  trees.push_back(treeNew);

  // original tree leaves
  Int_t           run   = 0;
  Int_t           event = 0;
  Int_t           lumi  = 0;
  Int_t           nvtx  = 0;
  Float_t         rho   = 0;
  Float_t         pu_weight = 0;
  Float_t         lumiWeight = 0;
  Float_t         perEveW   = 0;
  Int_t           numGenLevel = 0;
  Int_t           firedZcontrol = 0;
  vector<float>   *electron_eta = 0;
  vector<float>   *electron_pt = 0;
  vector<float>   *electron_r9 = 0;
  vector<bool>    *isTagMediumEle = 0;
  vector<bool>    *isTagTightEle = 0;
  vector<bool>    *isTagHltSafeEle = 0;
  vector<bool>    *electron_matchHLT = 0;
  vector<bool>    *electron_matchZHLT = 0;
  vector<bool>    *electron_matchMC  = 0;
  vector<float>   *gamma_pt  = 0;
  vector<float>   *gamma_eta = 0;
  vector<float>   *gamma_r9 = 0;
  vector<int>     *gamma_presel   = 0;
  vector<int>     *gamma_fullsel  = 0;
  vector<int>     *gamma_nm1chiso = 0;
  vector<int>     *gamma_nm1phiso = 0;
  vector<int>     *gamma_nm1hoe = 0;
  vector<int>     *gamma_nm1sieie = 0;
  vector<float>   *gamma_eleveto = 0;
  vector<bool>    *gamma_matchHLT  = 0;
  vector<bool>    *gamma_matchZHLT  = 0;
  vector<bool>    *gamma_matchMC  = 0;
  vector<float>   *gamma_sieie = 0;
  vector<float>   *gamma_hoe = 0;
  vector<float>   *gamma_chiso = 0;
  vector<float>   *gamma_phoiso = 0;
  vector<float>   *gamma_corrphoiso = 0;
  vector<float>   *invMass    = 0;
  vector<float>   *invMassRaw = 0;
  vector<int>     *eleIndex   = 0;
  vector<int>     *gammaIndex = 0;

  // List of branches - original tree
  TBranch        *b_run; 
  TBranch        *b_event;
  TBranch        *b_lumi;
  TBranch        *b_nvtx;
  TBranch        *b_rho;
  TBranch        *b_pu_weight;
  TBranch        *b_lumiWeight;
  TBranch        *b_perEveW;
  TBranch        *b_numGenLevel;
  TBranch        *b_firedZcontrol;
  TBranch        *b_electron_eta;   //!   
  TBranch        *b_electron_pt;   //!   
  TBranch        *b_electron_r9;   //!   
  TBranch        *b_isTagMediumEle;  
  TBranch        *b_isTagTightEle;
  TBranch        *b_isTagHltSafeEle;  
  TBranch        *b_electron_matchHLT; 
  TBranch        *b_electron_matchZHLT; 
  TBranch        *b_electron_matchMC;
  TBranch        *b_gamma_pt;   //!  
  TBranch        *b_gamma_eta;   //!   
  TBranch        *b_gamma_r9;   //!   
  TBranch        *b_gamma_presel; 
  TBranch        *b_gamma_fullsel; 
  TBranch        *b_gamma_nm1chiso;
  TBranch        *b_gamma_nm1phiso;
  TBranch        *b_gamma_nm1hoe;
  TBranch        *b_gamma_nm1sieie;
  TBranch        *b_gamma_eleveto; 
  TBranch        *b_gamma_matchHLT;
  TBranch        *b_gamma_matchZHLT;
  TBranch        *b_gamma_matchMC;
  TBranch        *b_gamma_sieie;
  TBranch        *b_gamma_hoe;
  TBranch        *b_gamma_chiso;
  TBranch        *b_gamma_phoiso;
  TBranch        *b_gamma_corrphoiso;
  TBranch        *b_invMass; 
  TBranch        *b_invMassRaw; 
  TBranch        *b_eleIndex;  
  TBranch        *b_gammaIndex;

  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("run", &run, &b_run);
  treeOrig->SetBranchAddress("event", &event, &b_event);
  treeOrig->SetBranchAddress("lumi", &lumi, &b_lumi);
  treeOrig->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  treeOrig->SetBranchAddress("rho", &rho, &b_rho);
  treeOrig->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  treeOrig->SetBranchAddress("lumiWeight", &lumiWeight, &b_lumiWeight);
  treeOrig->SetBranchAddress("perEveW", &perEveW, &b_perEveW);
  treeOrig->SetBranchAddress("numGenLevel", &numGenLevel, &b_numGenLevel);
  treeOrig->SetBranchAddress("firedZcontrol",&firedZcontrol,&b_firedZcontrol);
  treeOrig->SetBranchAddress("electron_eta", &electron_eta, &b_electron_eta);
  treeOrig->SetBranchAddress("electron_pt", &electron_pt, &b_electron_pt);
  treeOrig->SetBranchAddress("electron_r9", &electron_r9, &b_electron_r9);
  treeOrig->SetBranchAddress("isTagHltSafeEle", &isTagHltSafeEle, &b_isTagHltSafeEle);
  treeOrig->SetBranchAddress("isTagMediumEle", &isTagMediumEle, &b_isTagMediumEle);
  treeOrig->SetBranchAddress("isTagTightEle", &isTagTightEle, &b_isTagTightEle);
  treeOrig->SetBranchAddress("electron_matchHLT", &electron_matchHLT, &b_electron_matchHLT);
  treeOrig->SetBranchAddress("electron_matchZHLT", &electron_matchZHLT, &b_electron_matchZHLT);
  treeOrig->SetBranchAddress("electron_matchMC", &electron_matchMC, &b_electron_matchMC);
  treeOrig->SetBranchAddress("gamma_pt", &gamma_pt, &b_gamma_pt);
  treeOrig->SetBranchAddress("gamma_eta", &gamma_eta, &b_gamma_eta);
  treeOrig->SetBranchAddress("gamma_r9", &gamma_r9, &b_gamma_r9);
  treeOrig->SetBranchAddress("gamma_presel", &gamma_presel, &b_gamma_presel);
  treeOrig->SetBranchAddress("gamma_fullsel", &gamma_fullsel, &b_gamma_fullsel);
  treeOrig->SetBranchAddress("gamma_nm1chiso", &gamma_nm1chiso, &b_gamma_nm1chiso);
  treeOrig->SetBranchAddress("gamma_nm1phiso", &gamma_nm1phiso, &b_gamma_nm1phiso);
  treeOrig->SetBranchAddress("gamma_nm1hoe", &gamma_nm1hoe, &b_gamma_nm1hoe);
  treeOrig->SetBranchAddress("gamma_nm1sieie", &gamma_nm1sieie, &b_gamma_nm1sieie);
  treeOrig->SetBranchAddress("gamma_eleveto", &gamma_eleveto, &b_gamma_eleveto);
  treeOrig->SetBranchAddress("gamma_matchHLT", &gamma_matchHLT, &b_gamma_matchHLT);        
  treeOrig->SetBranchAddress("gamma_matchZHLT", &gamma_matchZHLT, &b_gamma_matchZHLT);        
  treeOrig->SetBranchAddress("gamma_matchMC", &gamma_matchMC, &b_gamma_matchMC);        
  treeOrig->SetBranchAddress("gamma_sieie", &gamma_sieie, &b_gamma_sieie);        
  treeOrig->SetBranchAddress("gamma_hoe", &gamma_hoe, &b_gamma_hoe);        
  treeOrig->SetBranchAddress("gamma_chiso", &gamma_chiso, &b_gamma_chiso);        
  treeOrig->SetBranchAddress("gamma_phoiso", &gamma_phoiso, &b_gamma_phoiso);        
  treeOrig->SetBranchAddress("gamma_corrphoiso", &gamma_corrphoiso, &b_gamma_corrphoiso);        
  treeOrig->SetBranchAddress("invMass", &invMass, &b_invMass);
  treeOrig->SetBranchAddress("invMassRaw", &invMassRaw, &b_invMassRaw);
  treeOrig->SetBranchAddress("eleIndex", &eleIndex, &b_eleIndex);
  treeOrig->SetBranchAddress("gammaIndex", &gammaIndex, &b_gammaIndex);

  // New variables
  float tag_pt, tag_r9;
  float tag_eta, tag_absEta;
  int tag_matchMC;
  int tag_matchHLT;
  int tag_matchZHLT;
  float probe_pt, probe_absEta;
  float probe_eta;
  float probe_r9;
  int probe_matchMC; 
  int probe_matchHLT; 
  int probe_matchZHLT; 
  int probe_fullsel;
  int probe_presel;
  int probe_tunedphisosel;
  int probe_ptonlycorrphisosel;
  int probe_nm1chiso;
  int probe_nm1phiso;
  int probe_nm1hoe;
  int probe_nm1sieie;
  int probe_totsel;
  int probe_eleveto;
  float probe_sieie, probe_hoe;
  float probe_chiso, probe_phoiso;
  float probe_corrphoiso;
  float mass;
  float massRaw;
  float xsecWeight, weight;
  
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];

    // New branches
    theTreeNew->Branch("run", &run, "run/I");
    theTreeNew->Branch("event", &event, "event/I");
    theTreeNew->Branch("lumi", &lumi, "lumi/I");
    theTreeNew->Branch("nvtx", &nvtx, "nvtx/I");
    theTreeNew->Branch("rho", &rho, "rho/F");
    theTreeNew->Branch("pu_weight", &pu_weight, "pu_weight/F");
    theTreeNew->Branch("firedZcontrol", &firedZcontrol, "firedZcontrol/I");
    theTreeNew->Branch("tag_pt",&tag_pt,"tag_pt/F");
    theTreeNew->Branch("tag_r9",&tag_r9,"tag_r9/F");
    theTreeNew->Branch("tag_eta",&tag_eta,"tag_eta/F");
    theTreeNew->Branch("tag_absEta",&tag_absEta,"tag_absEta/F");
    theTreeNew->Branch("tag_matchMC",&tag_matchMC,"tag_matchMC/I");
    theTreeNew->Branch("tag_matchHLT",&tag_matchHLT,"tag_matchHLT/I");
    theTreeNew->Branch("tag_matchZHLT",&tag_matchZHLT,"tag_matchZHLT/I");
    theTreeNew->Branch("probe_pt",&probe_pt,"probe_pt/F");
    theTreeNew->Branch("probe_absEta",&probe_absEta,"probe_absEta/F");
    theTreeNew->Branch("probe_eta",&probe_eta,"probe_eta/F");
    theTreeNew->Branch("probe_r9",&probe_r9,"probe_r9/F");
    theTreeNew->Branch("probe_fullsel", &probe_fullsel, "probe_fullsel/I");
    theTreeNew->Branch("probe_presel", &probe_presel, "probe_presel/I");
    theTreeNew->Branch("probe_tunedphisosel", &probe_tunedphisosel, "probe_tunedphisosel/I");
    theTreeNew->Branch("probe_ptonlycorrphisosel", &probe_ptonlycorrphisosel, "probe_ptonlycorrphisosel/I");
    theTreeNew->Branch("probe_nm1chiso", &probe_nm1chiso, "probe_nm1chiso/I");
    theTreeNew->Branch("probe_nm1phiso", &probe_nm1phiso, "probe_nm1phiso/I");
    theTreeNew->Branch("probe_nm1hoe", &probe_nm1hoe, "probe_nm1hoe/I");
    theTreeNew->Branch("probe_nm1sieie", &probe_nm1sieie, "probe_nm1sieie/I");
    theTreeNew->Branch("probe_totsel", &probe_totsel, "probe_totsel/I");
    theTreeNew->Branch("probe_matchHLT",&probe_matchHLT,"probe_matchHLT/I");
    theTreeNew->Branch("probe_matchZHLT",&probe_matchZHLT,"probe_matchZHLT/I");
    theTreeNew->Branch("probe_matchMC",&probe_matchMC,"probe_matchMC/I");
    theTreeNew->Branch("probe_eleveto",&probe_eleveto,"probe_eleveto/I");
    theTreeNew->Branch("probe_sieie",&probe_sieie,"probe_sieie/F");
    theTreeNew->Branch("probe_hoe",&probe_hoe,"probe_hoe/F");
    theTreeNew->Branch("probe_chiso",&probe_chiso,"probe_chiso/F");
    theTreeNew->Branch("probe_phoiso",&probe_phoiso,"probe_phoiso/F");
    theTreeNew->Branch("probe_corrphoiso",&probe_corrphoiso,"probe_corrphoiso/F");
    theTreeNew->Branch("mass", &mass, "mass/F");
    theTreeNew->Branch("massRaw", &massRaw, "massRaw/F");
    theTreeNew->Branch("xsecWeight", &xsecWeight, "xsecWeight/F");
    theTreeNew->Branch("weight", &weight, "weight/F");
  }

  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    for (unsigned int ii=0; ii<invMass->size(); ii++) {
      
      mass = (float)(invMass->at(ii));
      // if (mass<70 || mass>110) continue;     // to select the events in the good mass region *before* scale/smearings corrections
      
      // further selection on tag 
      float absTagEta = fabs(electron_eta->at(eleIndex->at(ii)));
      if (electron_pt->at(eleIndex->at(ii))<30)         continue;
      if (absTagEta>2.1)                                continue;   // tighter than on data, where the effective cut is around 2.18 
      if (!isTagHltSafeEle->at(eleIndex->at(ii)))       continue;
      if (!isTagTightEle->at(eleIndex->at(ii)))         continue;
      if (!electron_matchHLT->at(eleIndex->at(ii)))     continue;   
      
      // further selection on probe
      if (gamma_eleveto->at(gammaIndex->at(ii)))  continue;       // nominal T&P
      // if (!gamma_fullsel->at(gammaIndex->at(ii)))  continue;         // fake rate check  
      
      // only for HLT studies: the probe has to pass the full selection + inverted ele veto
      // if (!gamma_fullsel->at(gammaIndex->at(ii))) continue;       

      // match with mc-truth
      if (run==1){ 
	if (numGenLevel!=2 && numGenLevel!=0) cout << "Be careful: numGenLevel = " << numGenLevel << " - this is not contemplated!" << endl;
	if (numGenLevel==2 && !electron_matchMC->at(eleIndex->at(ii))) continue;   
	if (numGenLevel==2 && !gamma_matchMC->at(gammaIndex->at(ii)))  continue;
      }


      // To apply smearings and scales
      float tagEta   = electron_eta->at(eleIndex->at(ii));  
      float tagR9    = electron_r9->at(eleIndex->at(ii));
      float probeEta = gamma_eta->at(gammaIndex->at(ii));
      float probeR9  = gamma_r9->at(gammaIndex->at(ii)); 
      // smearings in MC
      if (run==1){
	float tagSmear   = smearings(tagR9, tagEta);
	float probeSmear = smearings(probeR9, probeEta);
	float theGaussSigma = 0.5 * sqrt (tagSmear*tagSmear + probeSmear*probeSmear);
	float theGaussMean  = 1.;
	float fromGauss = myRandom.Gaus(theGaussMean,theGaussSigma);
	mass = mass*fromGauss;
      }
      // scale corrections in data
      if (run>1){ 
	float tagScale = scaleCorr(run, tagR9, tagEta);
	float probeScale = scaleCorr(run, probeR9, probeEta);
	mass = sqrt(tagScale*probeScale)*mass;
      }

      // pass photon iso cut?
      float corrPhIso = gamma_corrphoiso->at(gammaIndex->at(ii));
      bool passPhIsoCut = isPhoIsoOk(probeEta, probeR9, corrPhIso);  

      // ad hoc cut to test photon isolation corrected for pT only, not for rho. 
      float probePt = gamma_pt->at(gammaIndex->at(ii));
      float ptOnlyCorrPhIso = gammaPtOnlyCorrPhIso(probeEta,probePt,gamma_phoiso->at(gammaIndex->at(ii)));
      bool passPtOnlyCorrPhIso = ptOnlyCorrPhIso<3.5;  // chiara: hardcoded. E' ragionevole per EE

      // invariant mass cut after smearings / scale correction
      if (mass<70 || mass>110) continue;

      // now making flat tree
      massRaw = (float)(invMassRaw->at(ii));
      tag_pt = electron_pt->at(eleIndex->at(ii));
      tag_r9 = electron_r9->at(eleIndex->at(ii));
      tag_eta = electron_eta->at(eleIndex->at(ii));
      tag_absEta = fabs(electron_eta->at(eleIndex->at(ii)));
      tag_matchHLT = electron_matchHLT->at(eleIndex->at(ii));
      tag_matchZHLT = electron_matchZHLT->at(eleIndex->at(ii));
      tag_matchMC = electron_matchMC->at(eleIndex->at(ii));
      probe_pt = gamma_pt->at(gammaIndex->at(ii));
      probe_absEta  = fabs(gamma_eta->at(gammaIndex->at(ii)));
      probe_eta  = gamma_eta->at(gammaIndex->at(ii));
      probe_r9  = gamma_r9->at(gammaIndex->at(ii));
      probe_fullsel = gamma_fullsel->at(gammaIndex->at(ii));  
      probe_presel = gamma_presel->at(gammaIndex->at(ii));  
      probe_tunedphisosel = passPhIsoCut && gamma_nm1phiso->at(gammaIndex->at(ii));                // full selection with tuned phiso cut
      probe_ptonlycorrphisosel = passPtOnlyCorrPhIso && gamma_nm1phiso->at(gammaIndex->at(ii));    // full selection replacing the ph-iso cut with an ad-=hoc one correcting for pt only
      probe_nm1chiso = gamma_nm1chiso->at(gammaIndex->at(ii));
      probe_nm1phiso = gamma_nm1phiso->at(gammaIndex->at(ii));
      probe_nm1hoe   = gamma_nm1hoe->at(gammaIndex->at(ii));
      probe_nm1sieie = gamma_nm1sieie->at(gammaIndex->at(ii));
      probe_totsel = gamma_fullsel->at(gammaIndex->at(ii)) && gamma_presel->at(gammaIndex->at(ii)); 
      probe_matchHLT = gamma_matchHLT->at(gammaIndex->at(ii));
      probe_matchZHLT = gamma_matchZHLT->at(gammaIndex->at(ii));
      probe_matchMC = gamma_matchMC->at(gammaIndex->at(ii));
      probe_eleveto = gamma_eleveto->at(gammaIndex->at(ii));
      probe_sieie = gamma_sieie->at(gammaIndex->at(ii)); 
      probe_hoe = gamma_hoe->at(gammaIndex->at(ii)); 
      probe_chiso = gamma_chiso->at(gammaIndex->at(ii)); 
      probe_phoiso = gamma_phoiso->at(gammaIndex->at(ii)); 
      probe_corrphoiso = gamma_corrphoiso->at(gammaIndex->at(ii)); 

      // weights
      if (run==1) {   // MC                                                                                                                   
	xsecWeight = perEveW * lumiForW * lumiWeight;
	weight     = xsecWeight * pu_weight;
      } else {
	xsecWeight = 1.;
	weight     = 1.;
      }

      treeNew->Fill();
      increm++;
    }
  }


  // new format
  cout << "treeNew = " << treeNew->GetEntries() << endl;
  cout << "increm = " << increm << endl;
  treeNew->Write();
  fileNew->Close();
  fileNew->ls();

  fileOrig->cd();
  fileOrig->Close();  
}

float scaleCorr(int run, float r9, float sceta) {
  
  float theScale = 1;

  if (fabs(sceta)<1 && r9>=0.94) {
    if (run>=273158 && run<=273301) theScale = 0.9952;
    if (run>=273302 && run<=273445) theScale = 0.9934;
    if (run>=273446 && run<=273491) theScale = 0.9937;
    if (run>=273492 && run<=273502) theScale = 0.9937;
    if (run>=273503 && run<=273727) theScale = 0.9933;
    if (run>=273728 && run<=274093) theScale = 0.9927;
    if (run>=274094 && run<=274171) theScale = 0.9911;
    if (run>=274172 && run<=274199) theScale = 0.9921;
    if (run>=274200 && run<=274243) theScale = 0.9916;
    if (run>=274244 && run<=274249) theScale = 0.9913;
    if (run==274250) theScale = 0.9919;
    if (run>=274251 && run<=274285) theScale = 0.9915;
    if (run>=274286 && run<=274315) theScale = 0.9930;
    if (run==274316) theScale = 0.9922;
    if (run>=274317 && run<=327335) theScale = 0.9923;
    if (run>=274336 && run<=274338) theScale = 0.9913;
    if (run>=274339 && run<=274344) theScale = 0.9914;
    if (run>=274345 && run<=274419) theScale = 0.9923;
    if (run>=274420 && run<=274439) theScale = 0.9917;
    if (run>=274440 && run<=274441) theScale = 0.9922;
    if (run>=274442 && run<=274953) theScale = 0.9919;
    if (run>=274954 && run<=274968) theScale = 0.9917;
    if (run==274969) theScale = 0.9910;
    if (run>=274970 && run<=274997) theScale = 0.9915;
    if (run==274998) theScale = 0.9918;
    if (run==274999) theScale = 0.9913;
    if (run>=275000 && run<=275058) theScale = 0.9913;
    if (run>=275059 && run<=275067) theScale = 0.9922;
    if (run>=275068 && run<=275072) theScale = 0.9919;
    if (run>=275073 && run<=275123) theScale = 0.9916;
    if (run>=275124 && run<=275125) theScale = 0.9914;
    if (run>=275282 && run<=275308) theScale = 0.9917;
    if (run>=275309 && run<=275310) theScale = 0.9911;
    if (run>=275311 && run<=275337) theScale = 0.9917;
    if (run>=275338 && run<=275369) theScale = 0.9921;
    if (run>=275370 && run<=275375) theScale = 0.9922;
    if (run>=275376 && run<=275657) theScale = 0.9916;
    if (run>=275658 && run<=275775) theScale = 0.9919;
    if (run>=275776 && run<=275831) theScale = 0.9909;
    if (run>=275832 && run<=275840) theScale = 0.9911;
    if (run>=275841 && run<=275885) theScale = 0.9905;
    if (run>=275886 && run<=275910) theScale = 0.9911;
    if (run>=275911 && run<=275917) theScale = 0.9921;
    if (run>=275918 && run<=276061) theScale = 0.9916;
    if (run>=276062) theScale = 0.9917;
  } 

  if (fabs(sceta)<1 && r9<0.94) {
    if (run>=273158 && run<=273301) theScale = 0.9948;
    if (run>=273302 && run<=273445) theScale = 0.9931;
    if (run>=273446 && run<=273491) theScale = 0.9933;
    if (run>=273492 && run<=273502) theScale = 0.9933;
    if (run>=273503 && run<=273727) theScale = 0.9929;
    if (run>=273728 && run<=274093) theScale = 0.9923;
    if (run>=274094 && run<=274171) theScale = 0.9907;
    if (run>=274172 && run<=274199) theScale = 0.9917;
    if (run>=274200 && run<=274243) theScale = 0.9912;
    if (run>=274244 && run<=274249) theScale = 0.9909;
    if (run==274250) theScale = 0.9916;
    if (run>=274251 && run<=274285) theScale = 0.9911;
    if (run>=274286 && run<=274315) theScale = 0.9926;
    if (run==274316) theScale = 0.9918;
    if (run>=274317 && run<=327335) theScale = 0.9919;
    if (run>=274336 && run<=274338) theScale = 0.9910;
    if (run>=274339 && run<=274344) theScale = 0.9910;
    if (run>=274345 && run<=274419) theScale = 0.9919;
    if (run>=274420 && run<=274439) theScale = 0.9913;
    if (run>=274440 && run<=274441) theScale = 0.9918;
    if (run>=274442 && run<=274953) theScale = 0.9915;
    if (run>=274954 && run<=274968) theScale = 0.9913;
    if (run==274969) theScale = 0.9906;
    if (run>=274970 && run<=274997) theScale = 0.9911;
    if (run==274998) theScale = 0.9914;
    if (run==274999) theScale = 0.9910;
    if (run>=275000 && run<=275058) theScale = 0.9909;
    if (run>=275059 && run<=275067) theScale = 0.9919;
    if (run>=275068 && run<=275072) theScale = 0.9915;
    if (run>=275073 && run<=275123) theScale = 0.9912;
    if (run>=275124 && run<=275125) theScale = 0.9910;
    if (run>=275282 && run<=275308) theScale = 0.9913;
    if (run>=275309 && run<=275310) theScale = 0.9907;
    if (run>=275311 && run<=275337) theScale = 0.9913;
    if (run>=275338 && run<=275369) theScale = 0.9917;
    if (run>=275370 && run<=275375) theScale = 0.9918;
    if (run>=275376 && run<=275657) theScale = 0.9912;
    if (run>=275658 && run<=275775) theScale = 0.9915;
    if (run>=275776 && run<=275831) theScale = 0.9905;
    if (run>=275832 && run<=275840) theScale = 0.9907;
    if (run>=275841 && run<=275885) theScale = 0.9901;
    if (run>=275886 && run<=275910) theScale = 0.9907;
    if (run>=275911 && run<=275917) theScale = 0.9917;
    if (run>=275918 && run<=276061) theScale = 0.9912;
    if (run>=276062) theScale = 0.9913;
  }

  if (fabs(sceta)<1.4442 && fabs(sceta)>=1 && r9>=0.94) {
    if (run>=273158 && run<=273301) theScale = 0.9897;
    if (run>=273302 && run<=273445) theScale = 0.9961;
    if (run>=273446 && run<=273491) theScale = 0.9869;
    if (run>=273492 && run<=273502) theScale = 0.9832;
    if (run>=273503 && run<=273727) theScale = 0.9857;
    if (run>=273728 && run<=274093) theScale = 0.9844;
    if (run>=274094 && run<=274171) theScale = 0.9850;
    if (run>=274172 && run<=274199) theScale = 0.9838;
    if (run>=274200 && run<=274243) theScale = 0.9843;
    if (run>=274244 && run<=274249) theScale = 0.9833;
    if (run==274250) theScale = 0.9879;
    if (run>=274251 && run<=274285) theScale = 0.9873;
    if (run>=274286 && run<=274315) theScale = 0.9862;
    if (run==274316) theScale = 0.9876;
    if (run>=274317 && run<=327335) theScale = 0.9871;
    if (run>=274336 && run<=274338) theScale = 0.9861;
    if (run>=274339 && run<=274344) theScale = 0.9844;
    if (run>=274345 && run<=274419) theScale = 0.9853;
    if (run>=274420 && run<=274439) theScale = 0.9856;
    if (run>=274440 && run<=274441) theScale = 0.9863;
    if (run>=274442 && run<=274953) theScale = 0.9851;
    if (run>=274954 && run<=274968) theScale = 0.9862;
    if (run==274969) theScale = 0.9863;
    if (run>=274970 && run<=274997) theScale = 0.9834;
    if (run==274998) theScale = 0.9876;
    if (run==274999) theScale = 0.9859;
    if (run>=275000 && run<=275058) theScale = 0.9841;
    if (run>=275059 && run<=275067) theScale = 0.9882;
    if (run>=275068 && run<=275072) theScale = 0.9855;
    if (run>=275073 && run<=275123) theScale = 0.9851;
    if (run>=275124 && run<=275125) theScale = 0.9847;
    if (run>=275282 && run<=275308) theScale = 0.9880;
    if (run>=275309 && run<=275310) theScale = 0.9845;
    if (run>=275311 && run<=275337) theScale = 0.9848;
    if (run>=275338 && run<=275369) theScale = 0.9868;
    if (run>=275370 && run<=275375) theScale = 0.9863;
    if (run>=275376 && run<=275657) theScale = 0.9857;
    if (run>=275658 && run<=275775) theScale = 0.9854;
    if (run>=275776 && run<=275831) theScale = 0.9849;
    if (run>=275832 && run<=275840) theScale = 0.9834;
    if (run>=275841 && run<=275885) theScale = 0.9832;
    if (run>=275886 && run<=275910) theScale = 0.9859;
    if (run>=275911 && run<=275917) theScale = 0.9863;
    if (run>=275918 && run<=276061) theScale = 0.9848;
    if (run>=276062) theScale = 0.9861;
  } 

  if (fabs(sceta)<1.4442 && fabs(sceta)>=1 && r9<0.94) {
    if (run>=273158 && run<=273301) theScale = 0.9934;
    if (run>=273302 && run<=273445) theScale = 0.9898;
    if (run>=273446 && run<=273491) theScale = 0.9906;
    if (run>=273492 && run<=273502) theScale = 0.9869;
    if (run>=273503 && run<=273727) theScale = 0.9894;
    if (run>=273728 && run<=274093) theScale = 0.9880;
    if (run>=274094 && run<=274171) theScale = 0.9887;
    if (run>=274172 && run<=274199) theScale = 0.9875;
    if (run>=274200 && run<=274243) theScale = 0.9879;
    if (run>=274244 && run<=274249) theScale = 0.9869;
    if (run==274250) theScale = 0.9907;
    if (run>=274251 && run<=274285) theScale = 0.9910;
    if (run>=274286 && run<=274315) theScale = 0.9898;
    if (run==274316) theScale = 0.9912;
    if (run>=274317 && run<=327335) theScale = 0.9908;
    if (run>=274336 && run<=274338) theScale = 0.9897;
    if (run>=274339 && run<=274344) theScale = 0.9880;
    if (run>=274345 && run<=274419) theScale = 0.9889;
    if (run>=274420 && run<=274439) theScale = 0.9893;
    if (run>=274440 && run<=274441) theScale = 0.9899;
    if (run>=274442 && run<=274953) theScale = 0.9888;
    if (run>=274954 && run<=274968) theScale = 0.9898;
    if (run==274969) theScale = 0.9899;
    if (run>=274970 && run<=274997) theScale = 0.9870;
    if (run==274998) theScale = 0.9913;
    if (run==274999) theScale = 0.9896;
    if (run>=275000 && run<=275058) theScale = 0.9877;
    if (run>=275059 && run<=275067) theScale = 0.9919;
    if (run>=275068 && run<=275072) theScale = 0.9891;
    if (run>=275073 && run<=275123) theScale = 0.9888;
    if (run>=275124 && run<=275125) theScale = 0.9883;
    if (run>=275282 && run<=275308) theScale = 0.9917;
    if (run>=275309 && run<=275310) theScale = 0.9882;
    if (run>=275311 && run<=275337) theScale = 0.9884;
    if (run>=275338 && run<=275369) theScale = 0.9904;
    if (run>=275370 && run<=275375) theScale = 0.9900;
    if (run>=275376 && run<=275657) theScale = 0.9893;
    if (run>=275658 && run<=275775) theScale = 0.9891;
    if (run>=275776 && run<=275831) theScale = 0.9885;
    if (run>=275832 && run<=275840) theScale = 0.9870;
    if (run>=275841 && run<=275885) theScale = 0.9869;
    if (run>=275886 && run<=275910) theScale = 0.9896;
    if (run>=275911 && run<=275917) theScale = 0.9899;
    if (run>=275918 && run<=276061) theScale = 0.9884;
    if (run>=276062) theScale = 0.9897;
  } 

  if (fabs(sceta)<2 && fabs(sceta)>=1.566 && r9>=0.94) {
    if (run>=273158 && run<=273301) theScale = 0.9913;
    if (run>=273302 && run<=273445) theScale = 0.9885;
    if (run>=273446 && run<=273491) theScale = 0.9881;
    if (run>=273492 && run<=273502) theScale = 0.9880;
    if (run>=273503 && run<=273727) theScale = 0.9880;
    if (run>=273728 && run<=274093) theScale = 0.9884;
    if (run>=274094 && run<=274171) theScale = 0.9927;
    if (run>=274172 && run<=274199) theScale = 0.9952;
    if (run>=274200 && run<=274243) theScale = 0.9963;
    if (run>=274244 && run<=274249) theScale = 0.9919;
    if (run==274250) theScale = 0.9907;
    if (run>=274251 && run<=274285) theScale = 0.9983;
    if (run>=274286 && run<=274315) theScale = 0.9979;
    if (run==274316) theScale = 0.9968;
    if (run>=274317 && run<=327335) theScale = 0.9984;
    if (run>=274336 && run<=274338) theScale = 0.9965;
    if (run>=274339 && run<=274344) theScale = 0.9949;
    if (run>=274345 && run<=274419) theScale = 0.9962;
    if (run>=274420 && run<=274439) theScale = 0.9972;
    if (run>=274440 && run<=274441) theScale = 0.9986;
    if (run>=274442 && run<=274953) theScale = 0.9942;
    if (run>=274954 && run<=274968) theScale = 0.9973;
    if (run==274969) theScale = 0.9974;
    if (run>=274970 && run<=274997) theScale = 0.9955;
    if (run==274998) theScale = 0.9976;
    if (run==274999) theScale = 0.9962;
    if (run>=275000 && run<=275058) theScale = 0.9936;
    if (run>=275059 && run<=275067) theScale = 0.9972;
    if (run>=275068 && run<=275072) theScale = 0.9973;
    if (run>=275073 && run<=275123) theScale = 0.9935;
    if (run>=275124 && run<=275125) theScale = 0.9972;
    if (run>=275282 && run<=275308) theScale = 0.9981;
    if (run>=275309 && run<=275310) theScale = 0.9938;
    if (run>=275311 && run<=275337) theScale = 0.9927;
    if (run>=275338 && run<=275369) theScale = 0.9969;
    if (run>=275370 && run<=275375) theScale = 0.9957;
    if (run>=275376 && run<=275657) theScale = 0.9932;
    if (run>=275658 && run<=275775) theScale = 0.9062;
    if (run>=275776 && run<=275831) theScale = 0.9946;
    if (run>=275832 && run<=275840) theScale = 0.9927;
    if (run>=275841 && run<=275885) theScale = 0.9913;
    if (run>=275886 && run<=275910) theScale = 0.9948;
    if (run>=275911 && run<=275917) theScale = 0.9957;
    if (run>=275918 && run<=276061) theScale = 0.9924;
    if (run>=276062) theScale = 0.9963;
  } 
  
  if (fabs(sceta)<2 && fabs(sceta)>=1.566 && r9<0.94) {
    if (run>=273158 && run<=273301) theScale = 0.9978;
    if (run>=273302 && run<=273445) theScale = 0.9851;
    if (run>=273446 && run<=273491) theScale = 0.9946;
    if (run>=273492 && run<=273502) theScale = 0.9945;
    if (run>=273503 && run<=273727) theScale = 0.9945;
    if (run>=273728 && run<=274093) theScale = 0.9950;
    if (run>=274094 && run<=274171) theScale = 0.9993;
    if (run>=274172 && run<=274199) theScale = 1.0018;
    if (run>=274200 && run<=274243) theScale = 1.0029;
    if (run>=274244 && run<=274249) theScale = 0.9984;
    if (run==274250) theScale = 1.0043;
    if (run>=274251 && run<=274285) theScale = 1.0048;
    if (run>=274286 && run<=274315) theScale = 1.0045;
    if (run==274316) theScale = 1.0034;
    if (run>=274317 && run<=327335) theScale = 1.0050;
    if (run>=274336 && run<=274338) theScale = 1.0030;
    if (run>=274339 && run<=274344) theScale = 1.0015;
    if (run>=274345 && run<=274419) theScale = 1.0027;
    if (run>=274420 && run<=274439) theScale = 1.0038;
    if (run>=274440 && run<=274441) theScale = 1.0052;
    if (run>=274442 && run<=274953) theScale = 1.0008;
    if (run>=274954 && run<=274968) theScale = 1.0038;
    if (run==274969) theScale = 1.0040;
    if (run>=274970 && run<=274997) theScale = 1.0020;
    if (run==274998) theScale = 1.0042;
    if (run==274999) theScale = 1.0027;
    if (run>=275000 && run<=275058) theScale = 1.0001;
    if (run>=275059 && run<=275067) theScale = 1.0038;
    if (run>=275068 && run<=275072) theScale = 1.0039;
    if (run>=275073 && run<=275123) theScale = 1.0000;
    if (run>=275124 && run<=275125) theScale = 1.0038;
    if (run>=275282 && run<=275308) theScale = 1.0047;
    if (run>=275309 && run<=275310) theScale = 1.0004;
    if (run>=275311 && run<=275337) theScale = 0.9993;
    if (run>=275338 && run<=275369) theScale = 1.0035;
    if (run>=275370 && run<=275375) theScale = 1.0023;
    if (run>=275376 && run<=275657) theScale = 0.9998;
    if (run>=275658 && run<=275775) theScale = 1.0027;
    if (run>=275776 && run<=275831) theScale = 1.0011;
    if (run>=275832 && run<=275840) theScale = 0.9992;
    if (run>=275841 && run<=275885) theScale = 0.9978;
    if (run>=275886 && run<=275910) theScale = 1.0014;
    if (run>=275911 && run<=275917) theScale = 1.0022;
    if (run>=275918 && run<=276061) theScale = 0.9989;
    if (run>=276062) theScale = 1.0029;
  } 

  if (fabs(sceta)<2.5 && fabs(sceta)>=2 && r9>=0.94) {
    if (run>=273158 && run<=273301) theScale = 0.9984;
    if (run>=273302 && run<=273445) theScale = 0.9861;
    if (run>=273446 && run<=273491) theScale = 0.9881;
    if (run>=273492 && run<=273502) theScale = 0.9879;
    if (run>=273503 && run<=273727) theScale = 0.9895;
    if (run>=273728 && run<=274093) theScale = 0.9880;
    if (run>=274094 && run<=274171) theScale = 0.9941;
    if (run>=274172 && run<=274199) theScale = 0.9955;
    if (run>=274200 && run<=274243) theScale = 0.9938;
    if (run>=274244 && run<=274249) theScale = 0.9949;
    if (run==274250) theScale = 0.9948;
    if (run>=274251 && run<=274285) theScale = 0.9916;
    if (run>=274286 && run<=274315) theScale = 0.9920;
    if (run==274316) theScale = 0.9927;
    if (run>=274317 && run<=327335) theScale = 0.9924;
    if (run>=274336 && run<=274338) theScale = 0.9912;
    if (run>=274339 && run<=274344) theScale = 0.9878;
    if (run>=274345 && run<=274419) theScale = 0.9911;
    if (run>=274420 && run<=274439) theScale = 0.9894;
    if (run>=274440 && run<=274441) theScale = 0.9920;
    if (run>=274442 && run<=274953) theScale = 0.9881;
    if (run>=274954 && run<=274968) theScale = 0.9911;
    if (run==274969) theScale = 0.9891;
    if (run>=274970 && run<=274997) theScale = 0.9893;
    if (run==274998) theScale = 0.9912;
    if (run==274999) theScale = 0.9880;
    if (run>=275000 && run<=275058) theScale = 0.9850;
    if (run>=275059 && run<=275067) theScale = 0.9926;
    if (run>=275068 && run<=275072) theScale = 0.9874;
    if (run>=275073 && run<=275123) theScale = 0.9859;
    if (run>=275124 && run<=275125) theScale = 0.9868;
    if (run>=275282 && run<=275308) theScale = 0.9890;
    if (run>=275309 && run<=275310) theScale = 0.9859;
    if (run>=275311 && run<=275337) theScale = 0.9844;
    if (run>=275338 && run<=275369) theScale = 0.9870;
    if (run>=275370 && run<=275375) theScale = 0.9868;
    if (run>=275376 && run<=275657) theScale = 0.9823;
    if (run>=275658 && run<=275775) theScale = 0.9901;
    if (run>=275776 && run<=275831) theScale = 0.9877;
    if (run>=275832 && run<=275840) theScale = 0.9855;
    if (run>=275841 && run<=275885) theScale = 0.9811;
    if (run>=275886 && run<=275910) theScale = 0.9840;
    if (run>=275911 && run<=275917) theScale = 0.9850;
    if (run>=275918 && run<=276061) theScale = 0.9848;
    if (run>=276062) theScale = 0.9840;
  } 

  if (fabs(sceta)<2.5 && fabs(sceta)>=2 && r9<0.94) {
    if (run>=273158 && run<=273301) theScale = 1.0050;
    if (run>=273302 && run<=273445) theScale = 0.9927;
    if (run>=273446 && run<=273491) theScale = 0.9947;
    if (run>=273492 && run<=273502) theScale = 0.9945;
    if (run>=273503 && run<=273727) theScale = 0.9961;
    if (run>=273728 && run<=274093) theScale = 0.9946;
    if (run>=274094 && run<=274171) theScale = 1.0007;
    if (run>=274172 && run<=274199) theScale = 1.0022;
    if (run>=274200 && run<=274243) theScale = 1.0004;
    if (run>=274244 && run<=274249) theScale = 1.0016;
    if (run==274250) theScale = 1.0014;
    if (run>=274251 && run<=274285) theScale = 0.9982;
    if (run>=274286 && run<=274315) theScale = 0.9986;
    if (run==274316) theScale = 0.9993;
    if (run>=274317 && run<=327335) theScale = 0.9990;
    if (run>=274336 && run<=274338) theScale = 0.9978;
    if (run>=274339 && run<=274344) theScale = 0.9944;
    if (run>=274345 && run<=274419) theScale = 0.9978;
    if (run>=274420 && run<=274439) theScale = 0.9961;
    if (run>=274440 && run<=274441) theScale = 0.9986;
    if (run>=274442 && run<=274953) theScale = 0.9947;
    if (run>=274954 && run<=274968) theScale = 0.9978;
    if (run==274969) theScale = 0.9957;
    if (run>=274970 && run<=274997) theScale = 0.9960;
    if (run==274998) theScale = 0.9979;
    if (run==274999) theScale = 0.9946;
    if (run>=275000 && run<=275058) theScale = 0.9916;
    if (run>=275059 && run<=275067) theScale = 0.9993;
    if (run>=275068 && run<=275072) theScale = 0.9940;
    if (run>=275073 && run<=275123) theScale = 0.9925;
    if (run>=275124 && run<=275125) theScale = 0.9934;
    if (run>=275282 && run<=275308) theScale = 0.9956;
    if (run>=275309 && run<=275310) theScale = 0.9925;
    if (run>=275311 && run<=275337) theScale = 0.9910;
    if (run>=275338 && run<=275369) theScale = 0.9936;
    if (run>=275370 && run<=275375) theScale = 0.9934;
    if (run>=275376 && run<=275657) theScale = 0.9889;
    if (run>=275658 && run<=275775) theScale = 0.9968;
    if (run>=275776 && run<=275831) theScale = 0.9944;
    if (run>=275832 && run<=275840) theScale = 0.9922;
    if (run>=275841 && run<=275885) theScale = 0.9877;
    if (run>=275886 && run<=275910) theScale = 0.9905;
    if (run>=275911 && run<=275917) theScale = 0.9916;
    if (run>=275918 && run<=276061) theScale = 0.9914;
    if (run>=276062) theScale = 0.9906;
  } 

  if (theScale==1) cout << "run = " << run << ": something fishi, uncorrected scale" << endl;

  return theScale;
}

float smearings(float r9, float sceta) {
  
  float theSmearings = 0.;

  if (fabs(sceta)<1 && r9>=0.94)     theSmearings = 0.0075;
  else if (fabs(sceta)<1 && r9<0.94) theSmearings = 0.0091;
  else if (fabs(sceta)<1.5 && fabs(sceta)>=1 && r9>=0.94) theSmearings = 0.0108;
  else if (fabs(sceta)<1.5 && fabs(sceta)>=1 && r9<0.94)  theSmearings = 0.0152;
  else if (fabs(sceta)<2 && fabs(sceta)>=1.5 && r9>=0.94) theSmearings = 0.0178;
  else if (fabs(sceta)<2 && fabs(sceta)>=1.5 && r9<0.94)  theSmearings = 0.0201;
  else if (fabs(sceta)<2.5 && fabs(sceta)>=2 && r9>=0.94) theSmearings = 0.0211;
  else if (fabs(sceta)<2.5 && fabs(sceta)>=2 && r9<0.94)  theSmearings = 0.0231;
  else cout << "something fishi, uncorrected resolution" << endl;

  return theSmearings;
}

bool isPhoIsoOk( float sceta, float r9, float corrphoiso) {

  // classes: 0 = EB highR9, 1 = EB low R9, 2 = EE high R9, 3 = EE lowR9                                                                                                                                       
  int etaclass = fabs(sceta)>1.5;
  int r9class  = r9<0.94;
  int theclass = 2.*etaclass + r9class;

  // cuts
  float phoiso_cut[4] = { 2.75, 2.75, 2., 2. };
  // float phoiso_cut[4] = { 3., 3., 2., 2. };
  // float phoiso_cut[4] = { 3.5, 3.5, 2., 2. };
  // float phoiso_cut[4] = { 4., 4., 2., 2. };

  if (corrphoiso > phoiso_cut[theclass]) return false;

  return true;
}

// photon isolation correction for pT only, not for rho
float gammaPtOnlyCorrPhIso( float sceta, float pt, float phoiso) {

  float ptCorrPhIso = 0;

  // alpha values - hardcoded 
  float phIsoAlpha[5] = { 2.5,2.5,2.5,2.5,2.5 };

  // kappa values - hardcoded 
  float phIsoKappa[5]= { 0.0045,0.0045,0.0045,0.003,0.003 };

  // EA corrections 
  int theEAregion  = effectiveAreaGammaRegion(sceta);

  // pt-only correction
  ptCorrPhIso = phIsoAlpha[theEAregion] + phoiso - phIsoKappa[theEAregion]*pt;

  return ptCorrPhIso;
}

int effectiveAreaGammaRegion(float sceta) {
  
  int theEAregion = 999;
  if (fabs(sceta)<=0.9) theEAregion = 0;
  if (fabs(sceta)<=1.5 && fabs(sceta)>0.9)  theEAregion = 1;
  if (fabs(sceta)<=2.0 && fabs(sceta)>1.5)  theEAregion = 2;
  if (fabs(sceta)<=2.2 && fabs(sceta)>2.0)  theEAregion = 3;
  if (fabs(sceta)<=2.5 && fabs(sceta)>2.2)  theEAregion = 4;
  return theEAregion;
}

