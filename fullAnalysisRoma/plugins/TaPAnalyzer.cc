#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"  

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"     
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"     
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"   
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"

#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

#define MAX_PU_REWEIGHT 60

using namespace std;
using namespace edm;
using namespace flashgg;
using namespace reco;

using pat::PackedGenParticle;   

class TaPAnalyzer : public edm::EDAnalyzer {
  
public:
  
  explicit TaPAnalyzer(const edm::ParameterSet&);
  ~TaPAnalyzer();
  
private:
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // weights
  void SetPuWeights(std::string puWeightFile);
  float GetPUWeight(float pun);
  
  // photons
  bool isGammaPresel( float sceta, float pt, float r9, float chiso);
  bool isGammaSelected( float rho, float pt, float sceta, float r9, float chiso, float nhiso, float phoiso, float hoe, float sieie, bool passElectronVeto);
  bool isGammaSelNm1ChIso( float rho, float pt, float sceta, float r9, float phoiso, float hoe, float sieie);
  bool isGammaSelNm1PhIso( float rho, float pt, float sceta, float r9, float chiso, float hoe, float sieie);
  bool isGammaSelNm1HoE( float rho, float pt, float sceta, float r9, float chiso, float phoiso, float sieie);
  bool isGammaSelNm1Sieie( float rho, float pt, float sceta, float r9, float chiso, float phoiso, float hoe );
  int effectiveAreaGammaRegion(float sceta);
  float gammaCorrPhIso( float rho, float pt, float sceta, float phoiso);
  
  // electrons
  float effectiveAreaEle03(float sceta);
  float effectiveArea2016Ele03(float sceta);
  Ptr<reco::Vertex> chooseElectronVertex( Ptr<flashgg::Electron> &elec, const std::vector<edm::Ptr<reco::Vertex> > &vertices );
  bool isMediumEle(float scEta, float hoe, float dphi, float deta, float sIeIe, float ep, float d0, float dz, float reliso, int missHits, bool passConvVeto) ;
  bool isTightEle(float scEta, float hoe, float dphi, float deta, float sIeIe, float ep, float d0, float dz, float reliso, int missHits, bool passConvVeto) ;
  bool isHLTSafeEle(float scEta, float hoe, float dphi, float detaOut, float sIeIe, float ep, float reliso, int missHits) ;

  void bookOutputTree();
  
  // collections
  EDGetTokenT<View<reco::Vertex> > vertexToken_;
  EDGetTokenT<View<Electron> > electronToken_;
  edm::EDGetTokenT<edm::View<flashgg::Photon> > photonToken_;
  EDGetTokenT<edm::View<PileupSummaryInfo> > PileUpToken_; 
  EDGetTokenT<double> rhoToken_; 
  EDGetTokenT<double> rhoEleToken_; 
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;  
  EDGetTokenT<GenEventInfoProduct> genInfoToken_;
  EDGetTokenT<View<pat::MET> > MetToken_;
  EDGetTokenT<View<reco::GenParticle> > genPartToken_;

  // sample-dependent parameters needed for the analysis
  int dopureweight_;
  string puWFileName_;
  int sampleIndex_;
  float lumiWeight_;
  
  // to compute weights for pileup
  std::vector<Double_t> puweights_;
  
  // to keep track of the number of events
  TH1F *h_entries;
  
  // events breakdown
  TH1F *h_selection;
  
  //---output tree branches variables
  edm::Service<TFileService> fs_;
  TTree* outTree_;
  
  int    run;
  int    event;
  int    lumi;
  int    nvtx;
  int    sampleID;
  float  rho;
  float  rhoEle;
  float  lumiWeight;
  float  pu_weight;
  float  pu_n;
  float perEveW;
  int numGenLevel;
  int firedZcontrol;
  int firedZtnp;

  float t1pfmet;
  
  int accEleSize;
  vector <float> electron_pt={};
  vector <float> electron_ene={};
  vector <float> electron_scRawEne={};
  vector <float> electron_eta={};
  vector <float> electron_phi={};
  vector <float> electron_r9={};
  vector <bool>  isTagHltSafeEle={};
  vector <bool>  isTagTightEle={};
  vector <bool>  isTagMediumEle={};
  vector <bool>  electron_matchHLT={};
  vector <bool>  electron_matchZHLT={};
  vector <bool>  electron_matchMC={};
  
  int  accGammaSize;                
  vector <float> gamma_pt={};
  vector <float> gamma_eta ={};
  vector <float> gamma_phi={};
  vector <float> gamma_r9={};
  vector <float> gamma_sieie={};
  vector <float> gamma_hoe ={};
  vector <float> gamma_scRawEne={};
  vector <float> gamma_ene={};
  vector <float> gamma_chiso={};
  vector <float> gamma_phoiso={};
  vector <float> gamma_corrphoiso={};
  vector <float> gamma_neuiso={};
  vector <float> gamma_eleveto={};
  vector <int>   gamma_presel={};
  vector <int>   gamma_fullsel={};
  vector <int>   gamma_nm1chiso={};
  vector <int>   gamma_nm1phiso={};
  vector <int>   gamma_nm1hoe={};
  vector <int>   gamma_nm1sieie={};
  vector <bool>  gamma_matchHLT={};
  vector <bool>  gamma_matchZHLT={};
  vector <bool>  gamma_matchMC={};
  
  vector <float> invMass={};
  vector <float> invMassRaw={};
  vector <int> eleIndex={};
  vector <int> gammaIndex={};
};  

TaPAnalyzer::TaPAnalyzer(const edm::ParameterSet& iConfig):

  // collections
  vertexToken_(consumes<View<reco::Vertex> >(iConfig.getUntrackedParameter<InputTag> ("VertexTag", InputTag("offlineSlimmedPrimaryVertices")))),
  electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag> ("ElectronTag"))),
  photonToken_(consumes<View<flashgg::Photon> >(iConfig.getUntrackedParameter<InputTag> ("PhotonTag", InputTag("flashggRandomizedPhotons")))),
  PileUpToken_(consumes<View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<InputTag> ("PileUpTag"))),
  rhoToken_(consumes<double>(iConfig.getParameter<InputTag> ("rhoTag"))), 
  rhoEleToken_(consumes<double>(iConfig.getParameter<InputTag> ("rhoEleTag"))), 
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  genInfoToken_( consumes<GenEventInfoProduct>(iConfig.getParameter<InputTag>("generatorInfo"))),
  MetToken_( consumes<View<pat::MET> >(iConfig.getParameter<InputTag>("MetTag" ))),
  genPartToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag> ("GenParticlesTag", InputTag("flashggPrunedGenParticles"))))
{ 
  dopureweight_ = iConfig.getUntrackedParameter<int>("dopureweight", 0);
  sampleIndex_  = iConfig.getUntrackedParameter<int>("sampleIndex",0);
  puWFileName_  = iConfig.getParameter<std::string>("puWFileName");   
  lumiWeight_   = iConfig.getUntrackedParameter<double>("lumiWeight",1.);
};

TaPAnalyzer::~TaPAnalyzer() { };

void TaPAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // To keep track of the total number of events
  h_entries->Fill(5);
  
  // Sample index
  sampleID   = sampleIndex_;
  lumiWeight = lumiWeight_;
  
  // --------------------------------------------------
  // access edm objects                                                                                    
  Handle<View<reco::Vertex> > primaryVertices;
  iEvent.getByToken(vertexToken_,primaryVertices);

  Handle<View<flashgg::Electron> > theElectrons;
  iEvent.getByToken( electronToken_, theElectrons );

  Handle<double> objs_rho;                                 
  iEvent.getByToken(rhoToken_,objs_rho);

  Handle<double> objs_rhoEle;     
  iEvent.getByToken(rhoEleToken_,objs_rhoEle);
  
  Handle<View<flashgg::Photon> > objs_pho;
  iEvent.getByToken(photonToken_,objs_pho);
  
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  Handle<View< PileupSummaryInfo> > PileupInfos;
  if (sampleID>0) iEvent.getByToken(PileUpToken_,PileupInfos);

  edm::Handle<GenEventInfoProduct> genInfo;      
  if (sampleID>0) iEvent.getByToken(genInfoToken_,genInfo);  

  edm::Handle<View<reco::GenParticle> > genParticles;
  if (sampleID>0) iEvent.getByToken( genPartToken_, genParticles );

  Handle<View<pat::MET> > METs;
  iEvent.getByToken( MetToken_, METs );
  if( METs->size() != 1 ) std::cout << "WARNING number of MET is not equal to 1" << std::endl; 
  Ptr<pat::MET> theMET = METs->ptrAt( 0 );
  t1pfmet = theMET->pt();
  

  // --------------------------------------------------
  // Event info
  run   = iEvent.eventAuxiliary().run();
  lumi  = iEvent.eventAuxiliary().luminosityBlock();
  event = iEvent.eventAuxiliary().event();

  // # Vertices
  nvtx = primaryVertices->size(); 
  
  // number of generated electrons
  numGenLevel = 0;

  // Energy density
  rho = *(objs_rho.product());
  rhoEle = *(objs_rhoEle.product());    
  
  // HLT paths preliminary settings
  firedZcontrol = 0;
  firedZtnp     = 0;

  // PU weight (for MC only and if requested)
  pu_weight = 1.;
  pu_n      = -1.;
  if (sampleID>0) {     // MC
    pu_n = 0.;
    for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI )
      {
	Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
	if( pu_bunchcrossing == 0 ) {
	  pu_n = PileupInfos->ptrAt( PVI )->getTrueNumInteractions();         
	}
      }
    if (dopureweight_) 
      pu_weight = GetPUWeight(pu_n);         
  }

  // other weights for the dataset
  perEveW = 1.;
  if (sampleID>0) { 
    const auto & eveWeights = genInfo->weights();
    if(!eveWeights.empty()) perEveW = eveWeights[0];
  }

  // Events breakdown  
  h_selection->Fill(0.,perEveW);
  

  // ----------------------------------------------------
  // save events only if:
  // 1) pass TnP HLT paths
  // 1) good vertex
  // 2) at least one tag
  // 3) at least one probe
  // ----------------------------------------------------
  
  
  // ----------------------------------------------------
  // 1) analysis cuts: trigger 

  // selected HLT object
  vector<float >hltTnpPt, hltTnpEta, hltTnpPhi;
  bool atLeastOneTag   = false;
  bool atLeastOneProbe = false;
    
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  
  // HLT paths for TnP
  string theTnPPath = "HLT_Ele27_eta2p1_WPLoose_Gsf_v";    // 2016 data: 8x MC is without HLT 

  // check if the event fired the TnP path
  bool fired = false;
  if (sampleID==0) {  // data only
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      string thisPath = names.triggerName(i);      
      //cout << "checkHLT: " << i << " " << thisPath << endl;
      if (thisPath.find(theTnPPath)==string::npos) continue;
      if (!triggerBits->accept(i)) continue;
      fired = true;
    }
  } else {
    fired = true;   // MC
  }
  
  if (sampleID==0) { 

    if (fired) {   
      
      h_selection->Fill(1.,perEveW);
      
      // HLT object firing the T&P path
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
	obj.unpackPathNames(names);
	
	vector<string> pathNamesAll = obj.pathNames(false);
	for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	  string thisPath = pathNamesAll[h];
	  
	  // the object has to be associated to the last filter of a succesfully path
	  bool isLF = obj.hasPathName( thisPath, true, false ); 
	  if (!isLF) continue;
	  
	  // the fired path must be our TnP path
	  if ( thisPath.find(theTnPPath)==string::npos) continue;
	  
	  hltTnpPt.push_back(obj.pt());
	  hltTnpEta.push_back(obj.eta());
	  hltTnpPhi.push_back(obj.phi());
	}	  
      }
    } 
  } else {

    h_selection->Fill(1.,perEveW);
  }

  if (sampleID==0) 
    if (hltTnpPt.size()!=hltTnpEta.size() || hltTnpPt.size()!=hltTnpPhi.size()) cout << "problem!" << endl;

  if ( hltTnpPt.size()>0 || sampleID>0 ) {  
    
    h_selection->Fill(2.,perEveW);
    
    // ----------------------------------------------------
    // 2) at least one good vertex found
    bool foundOne = false;
    for( unsigned int ivtx = 0 ; ivtx < primaryVertices->size() ; ivtx++ ) {
      if (foundOne) continue;
      
      bool goodVtx = true;
      float vtxX = primaryVertices->ptrAt(ivtx)->position().x();
      float vtxY = primaryVertices->ptrAt(ivtx)->position().y();
      float d0vtx = sqrt( vtxX*vtxX + vtxY*vtxY );
      if ( (primaryVertices->ptrAt(ivtx))->ndof()<=4 )  goodVtx = false;
      if ( fabs(d0vtx)>2 ) goodVtx = false;
      if ( fabs((primaryVertices->ptrAt(ivtx))->position().z())>=24 ) goodVtx = false;
      bool isVtxFake = ((primaryVertices->ptrAt(ivtx))->ndof()==0) && ((primaryVertices->ptrAt(ivtx))->chi2()==0);  
      if (isVtxFake) goodVtx = false;
      
      if (goodVtx) foundOne = true;
    } 
    if (foundOne) {
      
      h_selection->Fill(3.,perEveW);


      // Gen level match
      bool genEleFound = false;
      bool genPosFound = false;
      TLorentzVector myGenEle(0,0,0,0);  
      TLorentzVector myGenPos(0,0,0,0);  
      if (sampleID>0) {   
	for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
	  int status = genParticles->ptrAt( genLoop )->status();
	  int pdgid  = genParticles->ptrAt( genLoop )->pdgId();
	  if ( abs(pdgid)==11 && status==23 ) {
	    if (genParticles->ptrAt( genLoop )->mother(0)) {
	      if (genParticles->ptrAt( genLoop )->mother(0)->pdgId()==23) {
		float ptgen  = genParticles->ptrAt( genLoop )->pt();
		float etagen = genParticles->ptrAt( genLoop )->eta();
		float phigen = genParticles->ptrAt( genLoop )->phi();
		if (pdgid==11)  {
		  myGenPos.SetPtEtaPhiM(ptgen, etagen, phigen, 0.);
		  genPosFound = true;
		}
		if (pdgid==-11) {
		  myGenEle.SetPtEtaPhiM(ptgen, etagen, phigen, 0.);
		  genEleFound = true;
		}
		numGenLevel++;
	      }}}}
      }
    

      // ----------------------------------------------------  
      // Parenthesis: saving HLT infos for the Z control path
      string theZcontrolPath = "HLT_Ele35_WPLoose_Gsf_v";        
      vector<float >hltControlPt, hltControlEta, hltControlPhi;

      // check if the event fired the Z control pat
      firedZcontrol=0;
      if (sampleID==0) {  // data only
	for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
	  string thisPath = names.triggerName(i);      
	  if (thisPath.find(theZcontrolPath)==string::npos) continue;
	  if (!triggerBits->accept(i)) continue;
	  firedZcontrol = 1;
	}
      } else {      // MC
	firedZcontrol = 1;   
      }
      if (sampleID==0) {

	// HLT objects firing the Z control path
	for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
	  obj.unpackPathNames(names);
	  vector<string> pathNamesAll = obj.pathNames(false);
	  for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	    string thisPath = pathNamesAll[h];
	    bool isLF = obj.hasPathName( thisPath, true, false ); 
	    if (!isLF) continue;
	    if ( thisPath.find(theZcontrolPath)==string::npos) continue;
	    hltControlPt.push_back(obj.pt());
	    hltControlEta.push_back(obj.eta());
	    hltControlPhi.push_back(obj.phi());
	  }	  
	}
      }
    

      // ----------------------------------------------------  
      // 3) at least one tag candidate
      atLeastOneTag = false;
      
      // Loop over electron candidates
      const std::vector<edm::Ptr<reco::Vertex> > &vertexPointers = primaryVertices->ptrs();
      const std::vector<edm::Ptr<flashgg::Electron> > &ElectronPointers = theElectrons->ptrs();
      std::vector<int> acceptEle;
      // electrons in the acceptance
      for( unsigned int ElectronIndex = 0; ElectronIndex < ElectronPointers.size(); ElectronIndex++ ) {
	
	Ptr<flashgg::Electron> Electron = ElectronPointers[ElectronIndex];
	
	// acceptance
	float scEta = fabs( Electron->superCluster()->eta() );
	float elePt = Electron->pt();
	if( (fabs(scEta)>1.442 && fabs(scEta)<1.566) || fabs(scEta)>2.5 ) continue;
	if( elePt<20 ) continue;
	
	acceptEle.push_back(ElectronIndex);
      }
      
      // full cut based selection
      for (unsigned int iEle=0; iEle<acceptEle.size(); iEle++)
	{
	  int theOrigIndex = acceptEle.at(iEle);
	  Ptr<flashgg::Electron> Electron = ElectronPointers[theOrigIndex];
	  
	  // kine
	  float scEta  = fabs( Electron->superCluster()->eta() );
	  float elePt  = Electron->pt();
	  float eleEta = Electron->eta();
	  float elePhi = Electron->phi();

	  // this ele
	  TLorentzVector thisRecoEle(0,0,0,0);
	  thisRecoEle.SetPtEtaPhiM(elePt,eleEta,elePhi,0);

	  // match with selected HLT objects for the TnP path
	  bool matchHLT = false;
	  if (sampleID==0) {	    
	    for (int hltTnpC=0; hltTnpC<(int)hltTnpPt.size(); hltTnpC++)
	      {
		TLorentzVector thisHLTob(0,0,0,0);  
		float thisHLTpt  = hltTnpPt[hltTnpC];
		float thisHLTeta = hltTnpEta[hltTnpC];
		float thisHLTphi = hltTnpPhi[hltTnpC];
		thisHLTob.SetPtEtaPhiM(thisHLTpt,thisHLTeta,thisHLTphi,0);
		if(thisRecoEle.DeltaR(thisHLTob)<0.3)
		  matchHLT = true;
	      }
	  } else {
	    matchHLT = true;
	  }

	  // match with selected HLT objects for the Z control path
	  bool matchZHLT = false;
	  if (sampleID==0) {	    
	    for (int hltC=0; hltC<(int)hltControlPt.size(); hltC++)
	      {
		TLorentzVector thisHLTob(0,0,0,0);  
		float thisHLTpt  = hltControlPt[hltC];
		float thisHLTeta = hltControlEta[hltC];
		float thisHLTphi = hltControlPhi[hltC];
		thisHLTob.SetPtEtaPhiM(thisHLTpt,thisHLTeta,thisHLTphi,0);
		if(thisRecoEle.DeltaR(thisHLTob)<0.3)
		  matchZHLT = true;
	      }
	  } else {
	    matchZHLT = true;
	  }

	  // Match with MC truth
	  bool matchMC = false; 
	  if (sampleID>0) {  
	    if(genEleFound && thisRecoEle.DeltaR(myGenEle)<0.3) matchMC = true;  
	    if(genPosFound && thisRecoEle.DeltaR(myGenPos)<0.3) matchMC = true;  
	  }

	  // ID
	  float HoE = Electron->hcalOverEcal();
	  float DeltaPhiIn  = Electron->deltaPhiSuperClusterTrackAtVtx();
	  float DeltaEtaIn  = Electron->deltaEtaSuperClusterTrackAtVtx();
	  float DeltaEtaOut = Electron->deltaEtaSeedClusterTrackAtCalo();
	  float Full5x5Sieie = Electron->full5x5_sigmaIetaIeta(); 
	  float ecalEne = Electron->ecalEnergy();
	  float OneOverEoP;
	  if (ecalEne==0) {
	    cout << "electron energy is zero!! " << endl;
	    OneOverEoP = 1000000.;
	  } else {
	    OneOverEoP = 1.0/ecalEne - (Electron->eSuperClusterOverP())/ecalEne;
	  }
	  
	  // isolation with rho correction
	  reco::GsfElectron::PflowIsolationVariables pfIso = Electron->pfIsolationVariables();
	  // float corrHadPlusPho = pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rhoEle*effectiveArea2016Ele03(scEta);    
	  float corrHadPlusPho = pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rhoEle*effectiveAreaEle03(scEta);    
	  if (corrHadPlusPho<=0) corrHadPlusPho = 0.;
	  float absIsoWeffArea = pfIso.sumChargedHadronPt + corrHadPlusPho;
	  float relIso = absIsoWeffArea/elePt;
	  
	  // Impact parameters
	  Ptr<reco::Vertex> Electron_vtx = chooseElectronVertex( Electron, vertexPointers );
	  float d0 = Electron->gsfTrack()->dxy(Electron_vtx->position() );
	  float dz = Electron->gsfTrack()->dz( Electron_vtx->position() );
	  
	  // Conversions
	  const reco::HitPattern &hitPattern = Electron->gsfTrack()->hitPattern();
	  int mHits= hitPattern.numberOfHits(HitPattern::MISSING_INNER_HITS);
	  bool passConversionVeto = !(Electron->hasMatchedConversion());
	  
	  // All together
	  bool hltSafeEle = isHLTSafeEle(scEta, HoE, DeltaPhiIn, DeltaEtaOut, Full5x5Sieie, OneOverEoP, relIso, mHits) ;                             
	  bool mediumEle  = isMediumEle(scEta, HoE, DeltaPhiIn, DeltaEtaIn, Full5x5Sieie, OneOverEoP, d0, dz, relIso, mHits, passConversionVeto) ;
	  bool tightEle   = isTightEle(scEta, HoE, DeltaPhiIn, DeltaEtaIn, Full5x5Sieie, OneOverEoP, d0, dz, relIso, mHits, passConversionVeto) ;
	  
	  if (hltSafeEle && (mediumEle || tightEle)) atLeastOneTag = true;
	  
	  // Variables for the tree - for each electron in the acceptance - todo
	  electron_pt.push_back(elePt);
	  electron_ene.push_back(Electron->energy());
	  electron_scRawEne.push_back(Electron->superCluster()->rawEnergy());
	  electron_eta.push_back(Electron->superCluster()->eta());
	  electron_phi.push_back(Electron->superCluster()->phi());
	  electron_r9.push_back(Electron->full5x5_r9());    
	  isTagHltSafeEle.push_back(hltSafeEle);
	  isTagTightEle.push_back(tightEle);
	  isTagMediumEle.push_back(mediumEle);
	  electron_matchHLT.push_back(matchHLT);
	  electron_matchZHLT.push_back(matchZHLT);  
	  electron_matchMC.push_back(matchMC);                
	}  // tag
      accEleSize = electron_pt.size();
      
    
      // ----------------------------------------------------
      // 4) at least one probe found
      atLeastOneProbe = false;
      
      std::vector<int> acceptGamma;

      for(int phloop = 0; phloop < (int)objs_pho->size(); phloop++ ) {  
	
	Ptr<flashgg::Photon> g1 = objs_pho->ptrAt( phloop );
	
	// acceptance
	float gammaPt    = g1->et();
	float gammaScEta = (g1->superCluster())->eta();

	if (gammaPt<20) continue;
	if (fabs(gammaScEta)>2.5) continue;
	if (fabs(gammaScEta)>1.4442 && fabs(gammaScEta)<1.566) continue;
	
	acceptGamma.push_back(phloop);
      }

      // Photon candidates in the acceptance
      for(unsigned int iGamma=0; iGamma<acceptGamma.size(); iGamma++)
	{
	  int theOrigIndex = acceptGamma.at(iGamma);  
	  Ptr<flashgg::Photon> g1 = objs_pho->ptrAt( theOrigIndex );
	  
	  // kinematics
	  float scEta = (g1->superCluster())->eta();         
	  float pt    = g1->et();
	  float eta   = g1->eta();
	  float phi   = g1->phi();
	  TLorentzVector thisRecoGamma(0,0,0,0);
	  thisRecoGamma.SetPtEtaPhiM(pt,eta,phi,0);
	  
	  // match with selected HLT objects for the TnP path
	  bool matchHLT = false;
	  if (sampleID==0) {	    
	    for (int hltTnpC=0; hltTnpC<(int)hltTnpPt.size(); hltTnpC++)
	      {
		TLorentzVector thisHLTob(0,0,0,0);  
		float thisHLTpt  = hltTnpPt[hltTnpC];
		float thisHLTeta = hltTnpEta[hltTnpC];
		float thisHLTphi = hltTnpPhi[hltTnpC];
		thisHLTob.SetPtEtaPhiM(thisHLTpt,thisHLTeta,thisHLTphi,0);
		if(thisRecoGamma.DeltaR(thisHLTob)<0.3)
		  matchHLT = true;
	      }
	  } else {
	    matchHLT = true;
	  }

	  // match with selected HLT objects for the Z control path
	  bool matchZHLT = false;
	  if (sampleID==0) {	    
	    for (int hltC=0; hltC<(int)hltControlPt.size(); hltC++)
	      {
		TLorentzVector thisHLTob(0,0,0,0);  
		float thisHLTpt  = hltControlPt[hltC];
		float thisHLTeta = hltControlEta[hltC];
		float thisHLTphi = hltControlPhi[hltC];
		thisHLTob.SetPtEtaPhiM(thisHLTpt,thisHLTeta,thisHLTphi,0);
		if(thisRecoGamma.DeltaR(thisHLTob)<0.3)
		  matchZHLT = true;
	      }
	  } else {
	    matchZHLT = true;
	  }


	  // Match with MC truth
	  bool matchMC = false; 
	  if (sampleID>0) {  
	    if(genEleFound && thisRecoGamma.DeltaR(myGenEle)<0.3) matchMC = true;  
	    if(genPosFound && thisRecoGamma.DeltaR(myGenPos)<0.3) matchMC = true;  
	  }

	  // preselection
	  float R9noZS    = g1->full5x5_r9();    
	  float HoE       = g1->hadTowOverEm();
	  float sieienoZS = g1->full5x5_sigmaIetaIeta();   
	  
	  float chIso  = g1->egChargedHadronIso();
	  float neuIso = g1->egNeutralHadronIso();
	  float phoIso = g1->egPhotonIso();
	  
	  bool  eleVeto = g1->passElectronVeto();
	  
	  bool passPresel = isGammaPresel( scEta, pt, R9noZS, chIso); 

	  // full selection
	  bool passFullSelel = isGammaSelected( rho, pt, scEta, R9noZS, chIso, neuIso, phoIso, HoE, sieienoZS, eleVeto);

	  // corrected photon isolation
	  float corrPhIso = gammaCorrPhIso( rho, pt, scEta, phoIso ); 

	  // N-1 selections
	  bool passNm1ChIso = isGammaSelNm1ChIso( rho, pt, scEta, R9noZS, phoIso, HoE, sieienoZS);
	  bool passNm1PhIso = isGammaSelNm1PhIso( rho, pt, scEta, R9noZS, chIso, HoE, sieienoZS);
	  bool passNm1HoE   = isGammaSelNm1HoE( rho, pt, scEta, R9noZS, chIso, phoIso, sieienoZS);
	  bool passNm1Sieie = isGammaSelNm1Sieie( rho, pt, scEta, R9noZS, chIso, phoIso, HoE );

	  atLeastOneProbe = true;   // denominator = reco photons in acceptance with ET>20      
	
	  gamma_pt.push_back(pt);
	  gamma_eta.push_back(scEta);
	  gamma_phi.push_back(g1->superCluster()->phi());
	  gamma_r9.push_back(g1->full5x5_r9());
	  gamma_sieie.push_back(g1->full5x5_sigmaIetaIeta());
	  gamma_hoe.push_back(g1->hadTowOverEm());
	  gamma_scRawEne.push_back(g1->superCluster()->rawEnergy());
	  gamma_ene.push_back(g1->energy());
	  gamma_chiso.push_back(g1->egChargedHadronIso());
	  gamma_phoiso.push_back(g1->egPhotonIso());
	  gamma_corrphoiso.push_back(corrPhIso);
	  gamma_neuiso.push_back(g1->egNeutralHadronIso());
	  gamma_eleveto.push_back(g1->passElectronVeto());
	  gamma_presel.push_back(passPresel);
	  gamma_fullsel.push_back(passFullSelel);
	  gamma_nm1chiso.push_back(passNm1ChIso);
	  gamma_nm1phiso.push_back(passNm1PhIso);
	  gamma_nm1hoe.push_back(passNm1HoE);
	  gamma_nm1sieie.push_back(passNm1Sieie);
	  gamma_matchHLT.push_back(matchHLT);       
	  gamma_matchZHLT.push_back(matchZHLT);     
	  gamma_matchMC.push_back(matchMC);

	} // probe

      if (acceptGamma.size()>0 && !atLeastOneProbe) cout << "chiara: abbiamo un problema..." << endl;

    } // vertex
  } // HLT    
  accGammaSize = gamma_pt.size();   
      
  //---invariant mass and pt ratio
  for(int iGam=0; iGam<accGammaSize; ++iGam) {
    for(int iEle=0; iEle<accEleSize; ++iEle)  {
      
      TLorentzVector thisGamma(0,0,0,0); 
      TLorentzVector thisEle(0,0,0,0);  
      thisGamma.SetPtEtaPhiM(gamma_pt[iGam],gamma_eta[iGam],gamma_phi[iGam],0); 
      thisEle.SetPtEtaPhiM(electron_pt[iEle],electron_eta[iEle],electron_phi[iEle],0); 
      invMass.push_back((thisGamma+thisEle).M()); 

      TLorentzVector thisGammaRaw(0,0,0,0); 
      TLorentzVector thisEleRaw(0,0,0,0); 
      float rawPtGamma = gamma_pt[iGam]*gamma_scRawEne[iGam]/gamma_ene[iGam];
      thisGammaRaw.SetPtEtaPhiM(rawPtGamma,gamma_eta[iGam],gamma_phi[iGam],0); 
      float rawPtEle = electron_pt[iEle]*electron_scRawEne[iGam]/electron_ene[iGam];
      thisEleRaw.SetPtEtaPhiM(rawPtEle,electron_eta[iEle],electron_phi[iEle],0); 
      invMassRaw.push_back((thisGammaRaw+thisEleRaw).M()); 

      eleIndex.push_back(iEle);   
      gammaIndex.push_back(iGam);    
    }
  }
  
  //---fill output tree and reset
  if (atLeastOneProbe) {
    h_selection->Fill(4.,perEveW);
    if (atLeastOneTag) {
      h_selection->Fill(5.,perEveW);
      outTree_->Fill();
    }
  }
      
  //---tag
  electron_pt.clear();
  electron_ene.clear();
  electron_scRawEne.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_r9.clear();
  isTagHltSafeEle.clear();
  isTagTightEle.clear();
  isTagMediumEle.clear();
  electron_matchHLT.clear();
  electron_matchZHLT.clear();
  electron_matchMC.clear();
  //---probe
  gamma_pt.clear();
  gamma_eta.clear();
  gamma_phi.clear();
  gamma_r9.clear();
  gamma_sieie.clear();
  gamma_hoe.clear();
  gamma_scRawEne.clear();
  gamma_ene.clear();
  gamma_chiso.clear();
  gamma_phoiso.clear();
  gamma_corrphoiso.clear();
  gamma_neuiso.clear();
  gamma_eleveto.clear();
  gamma_presel.clear();
  gamma_fullsel.clear();
  gamma_nm1chiso.clear();
  gamma_nm1phiso.clear();
  gamma_nm1hoe.clear();
  gamma_nm1sieie.clear();
  gamma_matchHLT.clear();
  gamma_matchZHLT.clear();
  gamma_matchMC.clear();

  //---invariant mass
  invMass.clear();
  invMassRaw.clear();
  eleIndex.clear();
  gammaIndex.clear();
}

void TaPAnalyzer::beginJob()
{
    // loading weights for pileup if needed
    if (dopureweight_) 
        SetPuWeights(puWFileName_);
  
    // to keep track of the original number of events
    h_entries = fs_->make<TH1F>("h_entries", "h_entries", 10,  0., 10.);
    h_entries->Sumw2();

    // for the event breakdown  
    h_selection = fs_->make<TH1F>("h_selection", "h_selection", 6, -0.5, 5.5);
    h_selection->Sumw2();

    bookOutputTree();
}

void TaPAnalyzer::bookOutputTree() 
{
    outTree_ = fs_->make<TTree>("TaPtree", "TaPtree");
    
    cout << "Booking branches" << endl;
    
    outTree_->Branch("run", &run, "run/I");
    outTree_->Branch("event", &event, "event/I");
    outTree_->Branch("lumi", &lumi, "lumi/I");
    outTree_->Branch("nvtx", &nvtx, "nvtx/I");
    outTree_->Branch("rho", &rho, "rho/F");
    outTree_->Branch("rhoEle", &rhoEle, "rhoEle/F");
    outTree_->Branch("sampleID", &sampleID, "sampleID/I");
    outTree_->Branch("lumiWeight", &lumiWeight, "lumiWeight/F");
    outTree_->Branch("pu_weight", &pu_weight, "pu_weight/F");
    outTree_->Branch("pu_n", &pu_n, "pu_n/F");
    outTree_->Branch("perEveW", &perEveW, "perEveW/F");
    outTree_->Branch("numGenLevel", &numGenLevel, "numGenLevel/I");
    outTree_->Branch("firedZcontrol", &firedZcontrol, "firedZcontrol/I");
    outTree_->Branch("firedZtnp", &firedZtnp, "firedZtnp/I");

    outTree_->Branch("t1pfmet", &t1pfmet, "t1pfmet/F");

    outTree_->Branch("accEleSize", &accEleSize, "accEleSize/I");
    outTree_->Branch("electron_pt", "std::vector<float>", &electron_pt);
    outTree_->Branch("electron_ene", "std::vector<float>", &electron_ene);
    outTree_->Branch("electron_scRawEne", "std::vector<float>", &electron_scRawEne);
    outTree_->Branch("electron_eta", "std::vector<float>", &electron_eta);
    outTree_->Branch("electron_phi", "std::vector<float>", &electron_phi);
    outTree_->Branch("electron_r9", "std::vector<float>", &electron_r9);
    outTree_->Branch("isTagHltSafeEle", "std::vector<bool>", &isTagHltSafeEle );
    outTree_->Branch("isTagTightEle", "std::vector<bool>", &isTagTightEle );
    outTree_->Branch("isTagMediumEle", "std::vector<bool>", &isTagMediumEle );
    outTree_->Branch("electron_matchHLT", "std::vector<bool>", &electron_matchHLT );
    outTree_->Branch("electron_matchZHLT", "std::vector<bool>", &electron_matchZHLT );
    outTree_->Branch("electron_matchMC", "std::vector<bool>", &electron_matchMC );
 
    outTree_->Branch("accGammaSize",  &accGammaSize,  "accGammaSize/I");   
    outTree_->Branch("gamma_pt", "std::vector<float>", &gamma_pt);
    outTree_->Branch("gamma_eta", "std::vector<float>", &gamma_eta);
    outTree_->Branch("gamma_phi", "std::vector<float>", &gamma_phi);
    outTree_->Branch("gamma_r9", "std::vector<float>", &gamma_r9);
    outTree_->Branch("gamma_sieie", "std::vector<float>", &gamma_sieie);
    outTree_->Branch("gamma_hoe", "std::vector<float>", &gamma_hoe);
    outTree_->Branch("gamma_scRawEne", "std::vector<float>", &gamma_scRawEne);
    outTree_->Branch("gamma_ene", "std::vector<float>", &gamma_ene);
    outTree_->Branch("gamma_chiso", "std::vector<float>", &gamma_chiso);
    outTree_->Branch("gamma_phoiso", "std::vector<float>", &gamma_phoiso);
    outTree_->Branch("gamma_corrphoiso", "std::vector<float>", &gamma_corrphoiso);
    outTree_->Branch("gamma_neuiso", "std::vector<float>", &gamma_neuiso);
    outTree_->Branch("gamma_eleveto", "std::vector<float>", &gamma_eleveto);
    outTree_->Branch("gamma_presel", "std::vector<int>", &gamma_presel);
    outTree_->Branch("gamma_fullsel", "std::vector<int>", &gamma_fullsel);
    outTree_->Branch("gamma_nm1chiso", "std::vector<int>", &gamma_nm1chiso);
    outTree_->Branch("gamma_nm1phiso", "std::vector<int>", &gamma_nm1phiso);
    outTree_->Branch("gamma_nm1hoe", "std::vector<int>", &gamma_nm1hoe);
    outTree_->Branch("gamma_nm1sieie", "std::vector<int>", &gamma_nm1sieie);
    outTree_->Branch("gamma_matchHLT", "std::vector<bool>", &gamma_matchHLT );
    outTree_->Branch("gamma_matchZHLT", "std::vector<bool>", &gamma_matchZHLT );
    outTree_->Branch("gamma_matchMC", "std::vector<bool>", &gamma_matchMC );

    outTree_->Branch("invMass","std::vector<float>", &invMass);
    outTree_->Branch("invMassRaw","std::vector<float>", &invMassRaw);
    outTree_->Branch("eleIndex","std::vector<int>", &eleIndex);
    outTree_->Branch("gammaIndex","std::vector<int>", &gammaIndex);
}

void TaPAnalyzer::endJob() { }

// Phys14_25ns_v1
float TaPAnalyzer::effectiveAreaEle03(float theEta) {
  
    float theEA = -999;			   
    if(fabs(theEta) < 1) theEA = 0.1752; 
    else if(fabs(theEta) < 1.479) theEA = 0.1862;
    else if(fabs(theEta) < 2.0) theEA = 0.1411;
    else if(fabs(theEta) < 2.2) theEA = 0.1534;
    else if(fabs(theEta) < 2.3) theEA = 0.1903;
    else if(fabs(theEta) < 2.4) theEA = 0.2243;
    else if(fabs(theEta) < 2.5) theEA = 0.2687; 

    return theEA;
}

// Spring16 effective areas, see
// https://indico.cern.ch/event/482673/contributions/2187022/attachments/1282446/1905912/talk_electron_ID_spring16.pdf
float TaPAnalyzer::effectiveArea2016Ele03(float theEta) {
  
    float theEA = -999;			   
    if(fabs(theEta) < 1) theEA = 0.1703; 
    else if(fabs(theEta) < 1.479) theEA = 0.1715;
    else if(fabs(theEta) < 2.0) theEA = 0.1213;
    else if(fabs(theEta) < 2.2) theEA = 0.1230;
    else if(fabs(theEta) < 2.3) theEA = 0.1635;
    else if(fabs(theEta) < 2.4) theEA = 0.1937;
    else if(fabs(theEta) < 2.5) theEA = 0.2393; 

    return theEA;
}

// Electron vertex choice
Ptr<reco::Vertex> TaPAnalyzer::chooseElectronVertex( Ptr<flashgg::Electron> &elec, const std::vector<edm::Ptr<reco::Vertex> > &vertices ) {
  
    double vtx_dz = 1000000;
    unsigned int min_dz_vtx = -1;
    for( unsigned int vtxi = 0; vtxi < vertices.size(); vtxi++ ) {
        Ptr<reco::Vertex> vtx = vertices[vtxi];
        if( vtx_dz > fabs(elec->gsfTrack()->dz(vtx->position())) ) {
            vtx_dz = fabs( elec->gsfTrack()->dz( vtx->position() ) );
            min_dz_vtx = vtxi;
        }
    }
    return vertices[min_dz_vtx];
}

// Egamma Cut based ID 2015
bool TaPAnalyzer::isMediumEle(float scEta, float hoe, float dphi, float deta, float sIeIe,
			      float ep, float d0, float dz, float reliso, int missHits, bool passConvVeto) 
{
  
  bool okDeta, okDphi, okSieIe, okHoE, okEp, okD0, okDz, okIso, okMH, okConv;
  
  if (fabs(scEta)<1.479)
    { 
      okSieIe = sIeIe < 0.0101;
      okDeta  = fabs(deta) < 0.0103;
      okDphi  = fabs(dphi) < 0.0336;
      okHoE   = hoe < 0.0876;
      okIso   = reliso < 0.0766;
      okEp    = fabs(ep) < 0.0174;
      okD0    = fabs(d0) < 0.0118;
      okDz    = fabs(dz) < 0.373;
      okMH    = missHits<=2;
      okConv  = passConvVeto;
      
    } else {
    
    okSieIe = sIeIe < 0.0283;
    okDeta  = fabs(deta) < 0.00733;
    okDphi  = fabs(dphi) < 0.114;
    okHoE   = hoe < 0.0678;
    okIso   = reliso < 0.0678;
    okEp    = fabs(ep) < 0.0898;
    okD0    = fabs(d0) < 0.0739;
    okDz    = fabs(dz) < 0.602;
    okMH    = missHits<=1;
    okConv  = passConvVeto;
  }

  // chiara: for the moment we remove the deltaEta cut in EE
  bool okFullSel;
  // if (fabs(scEta)<1.479) 
  okFullSel = okDeta && okDphi && okSieIe && okHoE && okEp && okD0 && okDz && okIso && okMH && okConv;
  // else
  // okFullSel = okDphi && okSieIe && okHoE && okEp && okD0 && okDz && okIso && okMH && okConv;    

  return okFullSel;
}

// Egamma Cut based ID 2015
bool TaPAnalyzer::isTightEle(float scEta, float hoe, float dphi, float deta, float sIeIe,
			     float ep, float d0, float dz, float reliso, int missHits, bool passConvVeto)
{ 
 
    bool okDeta, okDphi, okSieIe, okHoE, okEp, okD0, okDz, okIso, okMH, okConv;
                                                                                                                                      
    if (fabs(scEta)<1.479)
      {
	okSieIe = sIeIe < 0.0101;
	okDeta  = fabs(deta) < 0.00926;
	okDphi  = fabs(dphi) < 0.0336;
	okHoE   = hoe < 0.0597;
	okIso   = reliso < 0.0354;
	okEp    = fabs(ep) < 0.012;
	okD0    = fabs(d0) < 0.0111;
	okDz    = fabs(dz) < 0.0466;
	okMH    = missHits<=2;
	okConv  = passConvVeto;
      }
    else
      {        
	okSieIe = sIeIe < 0.0279;
	okDeta  = fabs(deta) < 0.00724;
	okDphi  = fabs(dphi) < 0.0918;
	okHoE   = hoe < 0.0615;
	okIso   = reliso < 0.0646;
	okEp    = fabs(ep) < 0.00999;
	okD0    = fabs(d0) < 0.0351;
	okDz    = fabs(dz) < 0.417;
	okMH    = missHits<=1;
	okConv  = passConvVeto;
      }
    
    // chiara: for the moment we remove the deltaEta cut in EE
    bool okFullSel;
    // if (fabs(scEta)<1.479) 
    okFullSel = okDeta && okDphi && okSieIe && okHoE && okEp && okD0 && okDz && okIso && okMH && okConv;
    // else
    // okFullSel = okDphi && okSieIe && okHoE && okEp && okD0 && okDz && okIso && okMH && okConv;    

    return okFullSel;
}

// Egamma Cut based ID: 2016 HLT safe selection
// chiara: ci sara' una selezione on top di questa (retuning 2016) ma per ora uso questa + 2015 medium/tight
// https://indico.cern.ch/event/491536/contributions/2202104/attachments/1288921/1918602/talk_electron_ID_spring16.pdf
bool TaPAnalyzer::isHLTSafeEle(float scEta, float hoe, float dphi, float detaOut, float sIeIe, float ep, float reliso, int missHits) 
{  
  bool okDetaOut, okDphi, okSieIe, okHoE, okEp, okIso, okMH;
  
  if (fabs(scEta)<1.479)
    { 
      okSieIe   = sIeIe < 0.011;
      okDetaOut = fabs(detaOut) < 0.004;
      okDphi    = fabs(dphi) < 0.02;
      okEp      = fabs(ep) < 0.013;
      okIso     = reliso < 0.10;
      okHoE     = hoe < 0.06;
      okMH      = missHits<=1000000000;
      
    } else {
    
    okSieIe   = sIeIe < 0.031;
    okDetaOut = fabs(detaOut) < 0.007;
    okDphi    = fabs(dphi) < 0.020;
    okEp      = fabs(ep) < 0.013;
    okIso     = reliso < 0.10;
    okHoE     = hoe < 0.065;
    okMH      = missHits<=2;
  }
  
  bool okFullSel = okSieIe && okDetaOut && okDphi && okEp && okIso && okHoE && okMH;

  return okFullSel;
}

void TaPAnalyzer::SetPuWeights(std::string puWeightFile) {

    if (puWeightFile == "") {
        std::cout << "you need a weights file to use this function" << std::endl;
        return;
    }
    std::cout << "PU REWEIGHTING:: Using file " << puWeightFile << std::endl;

    TFile *f_pu  = new TFile(puWeightFile.c_str(),"READ");
    f_pu->cd();

    TH1D *puweights = 0;
    TH1D *gen_pu = 0;
    gen_pu    = (TH1D*) f_pu->Get("generated_pu");
    puweights = (TH1D*) f_pu->Get("weights");

    if (!puweights || !gen_pu) {
        std::cout << "weights histograms  not found in file " << puWeightFile << std::endl;
        return;
    }
    TH1D* weightedPU= (TH1D*)gen_pu->Clone("weightedPU");
    weightedPU->Multiply(puweights);

    // Rescaling weights in order to preserve same integral of events                               
    TH1D* weights = (TH1D*)puweights->Clone("rescaledWeights");
    weights->Scale( gen_pu->Integral(1,MAX_PU_REWEIGHT) / weightedPU->Integral(1,MAX_PU_REWEIGHT) );

    float sumPuWeights=0.;
    for (int i = 0; i<MAX_PU_REWEIGHT; i++) {
        float weight=1.;
        weight=weights->GetBinContent(i+1);
        sumPuWeights+=weight;
        puweights_.push_back(weight);
    }
}

float TaPAnalyzer::GetPUWeight(float pun) {
  
    float weight=1;
    if (sampleIndex_!=0 && pun<MAX_PU_REWEIGHT && puweights_.size()>0 && dopureweight_) 
        weight = puweights_[pun];
    return weight;
}

// miniAOD preselection + ECAL acceptance
bool TaPAnalyzer::isGammaPresel( float sceta, float pt, float r9, float chiso) {

    bool isPresel = false;

    // ECAL good acceptance
    if (fabs(sceta)>2.5) return false;
    if (fabs(sceta)>1.4442 && fabs(sceta)<1.566) return false;
  
    // miniAOD preselection
    if (r9>0.8)         return true;
    if (chiso<20)       return true;
    if ((chiso/pt)<0.3) return true;
  
    return isPresel;
}

bool TaPAnalyzer::isGammaSelected( float rho, float pt, float sceta, float r9, float chiso, float nhiso, float phoiso, float hoe, float sieie, bool passElectronVeto) {

    // classes: 0 = EB highR9, 1 = EB low R9, 2 = EE high R9, 3 = EE lowR9
    int etaclass = fabs(sceta)>1.5;
    int r9class  = r9<0.94;                   
    int theclass = 2.*etaclass + r9class;                  

    // cuts - hardcoded, v1
    float chiso_cut[4]     = { 5., 5., 5., 5. };
    float phoiso_cut[4]    = { 2.75, 2.75, 2., 2. };
    float sieie_cut[4]     = { 1.05e-02, 1.05e-02, 2.80e-02, 2.80e-02 };   
    float sieie_infCut[4]  = { 0.001, 0.001, 0.001, 0.001 };
    float hoe_cut[4]       = { 0.05, 0.05, 0.05, 0.05 };                   
  
    // effective areas - hardcoded 
    float phIsoAE[5] = { 0.17,0.14,0.11,0.14,0.22 };

    // alpha values - hardcoded
    float phIsoAlpha[5] = { 2.5,2.5,2.5,2.5,2.5 };

    // kappa values - hardcoded
    float phIsoKappa[5]= { 0.0045,0.0045,0.0045,0.003,0.003 };

    // EA corrections 
    int theEAregion  = effectiveAreaGammaRegion(sceta);   

    // full correction
    float corrPhIso = phIsoAlpha[theEAregion] + phoiso - rho*phIsoAE[theEAregion] - phIsoKappa[theEAregion]*pt;

    if (chiso > chiso_cut[theclass])      return false;
    if (corrPhIso > phoiso_cut[theclass]) return false;
    if (hoe > hoe_cut[theclass])          return false;
    if (sieie > sieie_cut[theclass])      return false;
    if (sieie < sieie_infCut[theclass])   return false;

    // electron veto 
    // if (!passElectronVeto) return false;

    return true;
} 

float TaPAnalyzer::gammaCorrPhIso( float rho, float pt, float sceta, float phoiso) {
  
  // effective areas - hardcoded 
  float phIsoAE[5] = { 0.17,0.14,0.11,0.14,0.22 };
  
  // alpha values - hardcoded
  float phIsoAlpha[5] = { 2.5,2.5,2.5,2.5,2.5 };
  
  // kappa values - hardcoded
  float phIsoKappa[5]= { 0.0045,0.0045,0.0045,0.003,0.003 };
  
  // EA corrections 
  int theEAregion  = effectiveAreaGammaRegion(sceta);   
  
  // full correction
  float corrPhIso = phIsoAlpha[theEAregion] + phoiso - rho*phIsoAE[theEAregion] - phIsoKappa[theEAregion]*pt;
  
  return corrPhIso;
} 

bool TaPAnalyzer::isGammaSelNm1ChIso( float rho, float pt, float sceta, float r9, float phoiso, float hoe, float sieie) {

    // classes: 0 = EB highR9, 1 = EB low R9, 2 = EE high R9, 3 = EE lowR9
    int etaclass = fabs(sceta)>1.5;
    int r9class  = r9<0.94;                   
    int theclass = 2.*etaclass + r9class;                  

    // cuts - hardcoded, v1
    float phoiso_cut[4]    = { 2.75, 2.75, 2., 2. };
    float sieie_cut[4]     = { 1.05e-02, 1.05e-02, 2.80e-02, 2.80e-02 };   
    float sieie_infCut[4]  = { 0.001, 0.001, 0.001, 0.001 };
    float hoe_cut[4]       = { 0.05, 0.05, 0.05, 0.05 };                   
  
    // effective areas - hardcoded 
    float phIsoAE[5] = { 0.17,0.14,0.11,0.14,0.22 };

    // alpha values - hardcoded
    float phIsoAlpha[5] = { 2.5,2.5,2.5,2.5,2.5 };

    // kappa values - hardcoded
    float phIsoKappa[5]= { 0.0045,0.0045,0.0045,0.003,0.003 };

    // EA corrections 
    int theEAregion  = effectiveAreaGammaRegion(sceta);   

    // full correction
    float corrPhIso = phIsoAlpha[theEAregion] + phoiso - rho*phIsoAE[theEAregion] - phIsoKappa[theEAregion]*pt;

    if (corrPhIso > phoiso_cut[theclass]) return false;
    if (hoe > hoe_cut[theclass])          return false;
    if (sieie > sieie_cut[theclass])      return false;
    if (sieie < sieie_infCut[theclass])   return false;

    return true;
} 

bool TaPAnalyzer::isGammaSelNm1PhIso( float rho, float pt, float sceta, float r9, float chiso, float hoe, float sieie) {

    // classes: 0 = EB highR9, 1 = EB low R9, 2 = EE high R9, 3 = EE lowR9
    int etaclass = fabs(sceta)>1.5;
    int r9class  = r9<0.94;                   
    int theclass = 2.*etaclass + r9class;                  

    // cuts - hardcoded, v1
    float chiso_cut[4]     = { 5., 5., 5., 5. };
    float sieie_cut[4]     = { 1.05e-02, 1.05e-02, 2.80e-02, 2.80e-02 };   
    float sieie_infCut[4]  = { 0.001, 0.001, 0.001, 0.001 };
    float hoe_cut[4]       = { 0.05, 0.05, 0.05, 0.05 };                   
  
    if (chiso > chiso_cut[theclass])      return false;
    if (hoe > hoe_cut[theclass])          return false;
    if (sieie > sieie_cut[theclass])      return false;
    if (sieie < sieie_infCut[theclass])   return false;

    return true;
} 

bool TaPAnalyzer::isGammaSelNm1HoE( float rho, float pt, float sceta, float r9, float chiso, float phoiso, float sieie) {

    // classes: 0 = EB highR9, 1 = EB low R9, 2 = EE high R9, 3 = EE lowR9
    int etaclass = fabs(sceta)>1.5;
    int r9class  = r9<0.94;                   
    int theclass = 2.*etaclass + r9class;                  

    // cuts - hardcoded, v1
    float chiso_cut[4]     = { 5., 5., 5., 5. };
    float phoiso_cut[4]    = { 2.75, 2.75, 2., 2. };
    float sieie_cut[4]     = { 1.05e-02, 1.05e-02, 2.80e-02, 2.80e-02 };   
    float sieie_infCut[4]  = { 0.001, 0.001, 0.001, 0.001 };
  
    // effective areas - hardcoded 
    float phIsoAE[5] = { 0.17,0.14,0.11,0.14,0.22 };

    // alpha values - hardcoded
    float phIsoAlpha[5] = { 2.5,2.5,2.5,2.5,2.5 };

    // kappa values - hardcoded
    float phIsoKappa[5]= { 0.0045,0.0045,0.0045,0.003,0.003 };

    // EA corrections 
    int theEAregion  = effectiveAreaGammaRegion(sceta);   

    // full correction
    float corrPhIso = phIsoAlpha[theEAregion] + phoiso - rho*phIsoAE[theEAregion] - phIsoKappa[theEAregion]*pt;

    if (chiso > chiso_cut[theclass])      return false;
    if (corrPhIso > phoiso_cut[theclass]) return false;
    if (sieie > sieie_cut[theclass])      return false;
    if (sieie < sieie_infCut[theclass])   return false;

    return true;
} 

bool TaPAnalyzer::isGammaSelNm1Sieie( float rho, float pt, float sceta, float r9, float chiso, float phoiso, float hoe ) {

    // classes: 0 = EB highR9, 1 = EB low R9, 2 = EE high R9, 3 = EE lowR9
    int etaclass = fabs(sceta)>1.5;
    int r9class  = r9<0.94;                   
    int theclass = 2.*etaclass + r9class;                  

    // cuts - hardcoded, v1
    float chiso_cut[4]     = { 5., 5., 5., 5. };
    float phoiso_cut[4]    = { 2.75, 2.75, 2., 2. };
    float hoe_cut[4]       = { 0.05, 0.05, 0.05, 0.05 };                   
  
    // effective areas - hardcoded 
    float phIsoAE[5] = { 0.17,0.14,0.11,0.14,0.22 };

    // alpha values - hardcoded
    float phIsoAlpha[5] = { 2.5,2.5,2.5,2.5,2.5 };

    // kappa values - hardcoded
    float phIsoKappa[5]= { 0.0045,0.0045,0.0045,0.003,0.003 };

    // EA corrections 
    int theEAregion  = effectiveAreaGammaRegion(sceta);   

    // full correction
    float corrPhIso = phIsoAlpha[theEAregion] + phoiso - rho*phIsoAE[theEAregion] - phIsoKappa[theEAregion]*pt;

    if (chiso > chiso_cut[theclass])      return false;
    if (corrPhIso > phoiso_cut[theclass]) return false;
    if (hoe > hoe_cut[theclass])          return false;

    return true;
} 

int TaPAnalyzer::effectiveAreaGammaRegion(float sceta) {

    int theEAregion = 999;
    if (fabs(sceta)<=0.9) theEAregion = 0;
    if (fabs(sceta)<=1.5 && fabs(sceta)>0.9)  theEAregion = 1;
    if (fabs(sceta)<=2.0 && fabs(sceta)>1.5)  theEAregion = 2;   
    if (fabs(sceta)<=2.2 && fabs(sceta)>2.0)  theEAregion = 3;
    if (fabs(sceta)<=2.5 && fabs(sceta)>2.2)  theEAregion = 4;
    return theEAregion;
}


DEFINE_FWK_MODULE(TaPAnalyzer);

