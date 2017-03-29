// -*- C++ -*-
//
// Package:    CEPDiJetsAnalyzer/AnalyzerPPS
// Class:      AnalyzerPPS
// 
/**\class AnalyzerPPS AnalyzerPPS.cc CEPDiJetsAnalyzer/AnalyzerPPS/plugins/AnalyzerPPS.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Pacheco
//         Created:  Tue, 07 Jun 2016 16:25:18 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Event Info
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// Gen Info
  // Gen Jets
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
  // Gen Particles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// Math
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

// Central Detector (CMS)
  // Jets
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
  // Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

  // Tracks Associated with Jets
#include "DataFormats/JetReco/interface/JetTrackMatch.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"

// CT-PPS
#include "FastSimulation/PPSFastObjects/interface/PPSSpectrometer.h"
#include "FastSimulation/PPSFastObjects/interface/PPSGenData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSSimData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSRecoData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSGenVertex.h"
#include "FastSimulation/PPSFastObjects/interface/PPSRecoVertex.h"

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// C/C++
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <map>
#include <iostream>


using namespace edm;
using namespace std;


//
// class declaration
//

class AnalyzerPPS : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AnalyzerPPS(const edm::ParameterSet&);
      ~AnalyzerPPS();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      void Init();
      void MCGenInfo               (const edm::Event&, const edm::EventSetup&);
      void GenCollections          (const edm::Event&, const edm::EventSetup&);
      void GenPPSInfo              (const edm::Event&, const edm::EventSetup&);
      void FillCollections         (const edm::Event&, const edm::EventSetup&);
      void SortingObjects          (const edm::Event&, const edm::EventSetup&);
      void AssociateJetsWithVertex (const edm::Event&, const edm::EventSetup&);



      
      typedef edm::View<reco::Vertex> VertexCollection;
      typedef PPSGenData    PPSGen; 
      typedef PPSSimData    PPSSim; 
      typedef PPSRecoData   PPSReco; 

      TTree* tree;

      double EBeam_;
//      double PPSres_;
//      double Nsigma_;

      // Boolean Flags
      bool runWithWeightGen_;

      // MCGenInfo function
      double GeneratorWeight, Pthat;
      int eventNumber, runNumber;

      // PPSGenInfo function
      double PPS_Gen_xiARMF, PPS_Gen_tARMF, PPS_Gen_etaARMF, PPS_Gen_phiARMF;
      double PPS_Gen_xiARMB, PPS_Gen_tARMB, PPS_Gen_etaARMB, PPS_Gen_phiARMB;
      std::vector<double> PPSGenVertexVectorX, PPSGenVertexVectorY, PPSGenVertexVectorZ;
     
      // GenCollection function
        // Gen Proton Info
      std::vector<const reco::GenParticle*> GenProtonVectorInfo;
      std::vector< math::XYZTLorentzVector > protonLorentzVector;
       //
      double GenMxx,GenMx, GenXi0, GenXi1, GenMjj, GenRjj;
      int Genjetsize;
      std::vector<const reco::GenJet*> GenJetsVector;
      std::vector<double> GenJetsVector_pt, GenJetsVector_eta, GenJetsVector_phi;

      // FillCollection function
	// vertex
      std::vector<const reco::Vertex*> VertexVector;
      std::vector<double> VertexVectorX, VertexVectorY, VertexVectorZ;
      int nTracks, nVertex;
	// jets
      std::vector<const pat::Jet*> JetsVector;
      double Mx;
	//dijets 
      TLorentzVector wj1, wj2, wdijet;
      double MJJWide, DeltaEtaJJWide, DeltaPhiJJWide, RJJWide, wj1Pt, wj2Pt;
	// pps reco
      std::vector<double> PPSRecoVertexVectorX, PPSRecoVertexVectorY, PPSRecoVertexVectorZ; 
      int PPSRecoVertexSize;
      std::vector< std::pair<double,double> > PPSCMSVertex; //unused
      int indexGold, nGV, nGV_sel;
      double diffVertex;
	// arm F
      std::vector<double> xiPPSARMF, tPPSARMF, etaPPSARMF, phiPPSARMF;
      std::vector<double> thetaXPPSARMF, thetaYPPSARMF ;
      std::vector<double> xDet1PPSARMF, yDet1PPSARMF, xDet2PPSARMF, yDet2PPSARMF;
      std::vector<double> tofPPSARMF, xTofPPSARMF, yTofPPSARMF;
	// arm B
      std::vector<double> xiPPSARMB, tPPSARMB, etaPPSARMB, phiPPSARMB;
      std::vector<double> thetaXPPSARMB, thetaYPPSARMB ;
      std::vector<double> xDet1PPSARMB, yDet1PPSARMB, xDet2PPSARMB, yDet2PPSARMB;
      std::vector<double> tofPPSARMB, xTofPPSARMB, yTofPPSARMB;
     
      // SortingObjects function
      std::vector<double> JetsVector_pt, JetsVector_eta, JetsVector_phi;
      double Mjj, Rjj;
      double yRapidity;
      int jetsize;
      double DeltaPhi12, DeltaEta12;

      // vertex associating function
      double GoldenVertexZ;
      std::vector<double> MinimumDistance;
      double MinDistance, MaxDistance;




      // Tokens
	// gen particles
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken;
	// gen jets
      edm::EDGetTokenT<reco::GenJetCollection> genjetsToken;
	// vertex
      edm::EDGetTokenT<VertexCollection> vertexToken;
	// tracks
//      edm::EDGetTokenT<TracksCollection> trackToken;
	// jets
      edm::EDGetTokenT<pat::JetCollection> jetsToken;
        // widejets
      edm::EDGetTokenT<std::vector<math::PtEtaPhiMLorentzVector> > widejetsToken;
	// pps gen
      edm::EDGetTokenT<PPSSpectrometer<PPSGen> > ppsGenToken;
	// pps sim
      edm::Handle<PPSSpectrometer<PPSSim> > ppsSim;
      edm::EDGetTokenT<PPSSpectrometer<PPSSim> > ppsSimToken;
	// pps reco
      edm::EDGetTokenT<PPSSpectrometer<PPSReco> > ppsRecoToken;

      
      
      // histogram
      TH1F *h_eventSelection;
	// Primary Vertex Information
      double pv_multiplicity;
      TH1F *h_pv_multiplicity , *h_pv_x	, *h_pv_y,  *h_pv_z;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AnalyzerPPS::AnalyzerPPS(const edm::ParameterSet& iConfig):
          genParticleToken (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag> ("genParticle")))
        , genjetsToken (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag> ("genjets")))
	, vertexToken (consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertex")))
	, jetsToken (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
	, widejetsToken (consumes<std::vector<math::PtEtaPhiMLorentzVector>>(iConfig.getUntrackedParameter<edm::InputTag>("widejetsTag")))
	, ppsGenToken (consumes<PPSSpectrometer<PPSGen>>(iConfig.getParameter<edm::InputTag>("ppsGen")))
//	, ppsSimToken (consumes<PPSSpectrometer<PPSSim>>(iConfig.getParameter<edm::InputTag>("ppsSim")))
	, ppsRecoToken (consumes<PPSSpectrometer<PPSReco>>(iConfig.getParameter<edm::InputTag>("ppsReco")))
{
   //now do what ever initialization is needed
   //usesResource("TFileService");

   edm::Service<TFileService> fs;

   runWithWeightGen_     = (iConfig.getParameter<bool>("RunWithWeightGen"));
   EBeam_                = (iConfig.getParameter<double>("EBeam"));
//   PPSres_               = (iConfig.getParameter<double>("PPSres"));;
//   Nsigma_               = (iConfig.getParameter<double>("Nsigma"));;

   tree = fs->make<TTree>( "PPS","PPS" );
   tree->Branch("pv_multiplicity",&pv_multiplicity,"pv_multiplicity/D");

   tree->Branch("GeneratorWeight",&GeneratorWeight,"GeneratorWeight/D"); 
   tree->Branch("ProtonsP4",&protonLorentzVector);
   tree->Branch("GenMx",&GenMx,"GenMx/D");
   tree->Branch("GenXi0",&GenXi0,"GenXi0/D"); 
   tree->Branch("GenXi1",&GenXi1,"GenXi1/D");   
   tree->Branch("GenJetsMultiplicity",&Genjetsize, "GenJetsMultiplicity/I");
   tree->Branch("GenJetsPt",&GenJetsVector_pt);
   tree->Branch("GenJetsEta",&GenJetsVector_eta);
   tree->Branch("GenJetsPhi",&GenJetsVector_phi);
   tree->Branch("GenMjj",&GenMjj,"GenMjj/D");
   tree->Branch("GenRjj",&GenRjj,"GenRjj/D");
   tree->Branch("PPS_Gen_xiARMF", &PPS_Gen_xiARMF, "PPS_Gen_xiARMF/D");
   tree->Branch("PPS_Gen_tARMF",  &PPS_Gen_tARMF, "PPS_Gen_tARMF/D");
   tree->Branch("PPS_Gen_etaARMF",&PPS_Gen_etaARMF, "PPS_Gen_etaARMF/D");
   tree->Branch("PPS_Gen_phiARMF",&PPS_Gen_phiARMF, "PPS_Gen_phiARMF/D");
   tree->Branch("PPS_Gen_xiARMB", &PPS_Gen_xiARMB, "PPS_Gen_xiARMB/D");
   tree->Branch("PPS_Gen_tARMB",  &PPS_Gen_tARMB, "PPS_Gen_tARMB/D");
   tree->Branch("PPS_Gen_etaARMB",&PPS_Gen_etaARMB, "PPS_Gen_etaARMB/D");
   tree->Branch("PPS_Gen_phiARMB",&PPS_Gen_phiARMB, "PPS_Gen_phiARMB/D");
   tree->Branch("PPSGenVertexVector_x",&PPSGenVertexVectorX);
   tree->Branch("PPSGenVertexVector_y",&PPSGenVertexVectorY);
   tree->Branch("PPSGenVertexVector_z",&PPSGenVertexVectorZ);
   tree->Branch("nVertex",&nVertex,"nVertex/I");
   tree->Branch("PPSRecoVertexSize",&PPSRecoVertexSize, "PPSRecoVertexSize/I");   
   tree->Branch("PPSRecoVertexVector_x",&PPSRecoVertexVectorX);
   tree->Branch("PPSRecoVertexVector_y",&PPSRecoVertexVectorY);
   tree->Branch("PPSRecoVertexVector_z",&PPSRecoVertexVectorZ);
   tree->Branch("xiPPSARMF", &xiPPSARMF);
   tree->Branch("tPPSARMF", &tPPSARMF);
   tree->Branch("etaPPSARMF", &etaPPSARMF);
   tree->Branch("phiPPSARMF", &phiPPSARMF);
   tree->Branch("thetaXPPSARMF", &thetaXPPSARMF);
   tree->Branch("thetaYPPSARMF", &thetaYPPSARMF);
   tree->Branch("tofPPSARMF", &tofPPSARMF);
   tree->Branch("xDet1PPSARMF", &xDet1PPSARMF);
   tree->Branch("yDet1PPSARMF", &yDet1PPSARMF);
   tree->Branch("xDet2PPSARMF", &xDet2PPSARMF);
   tree->Branch("yDet2PPSARMF", &yDet2PPSARMF);
   tree->Branch("xTofPPSARMF", &xTofPPSARMF);
   tree->Branch("yTofPPSARMF", &yTofPPSARMF);
   tree->Branch("xiPPSARMB", &xiPPSARMB);
   tree->Branch("tPPSARMB", &tPPSARMB);
   tree->Branch("etaPPSARMB", &etaPPSARMB);
   tree->Branch("phiPPSARMB", &phiPPSARMB);
   tree->Branch("thetaXPPSARMB", &thetaXPPSARMB);
   tree->Branch("thetaYPPSARMB", &thetaYPPSARMB);
   tree->Branch("tofPPSARMB", &tofPPSARMB);
   tree->Branch("xDet1PPSARMB", &xDet1PPSARMB);
   tree->Branch("yDet1PPSARMB", &yDet1PPSARMB);
   tree->Branch("xDet2PPSARMB", &xDet2PPSARMB);
   tree->Branch("yDet2PPSARMB", &yDet2PPSARMB);
   tree->Branch("xTofPPSARMB", &xTofPPSARMB);
   tree->Branch("yTofPPSARMB", &yTofPPSARMB);
   tree->Branch("Mx",&Mx,"Mx/D");
   tree->Branch("VertexVectorX",&VertexVectorX);
   tree->Branch("VertexVectorY",&VertexVectorY);
   tree->Branch("VertexVectorZ",&VertexVectorZ);
   tree->Branch("JetsMultiplicity", &jetsize, "JetsMultiplicity/I");
   tree->Branch("JetsPt",&JetsVector_pt);
   tree->Branch("JetsEta",&JetsVector_eta);
   tree->Branch("JetsPhi",&JetsVector_phi);
   tree->Branch("DeltaPhi12",&DeltaPhi12);
   tree->Branch("DeltaEta12",&DeltaEta12);
   tree->Branch("Mjj",&Mjj,"Mjj/D");
   tree->Branch("Rjj",&Rjj,"Rjj/D");
   tree->Branch("yRapidity",&yRapidity,"yRapidity/D");
   tree->Branch("GoldenVertexZ",&GoldenVertexZ,"GoldenVertexZ/D");
   tree->Branch("MinDistance",&MinDistance,"MinDistance/D");
   tree->Branch("MaxDistance",&MaxDistance,"MaxDistance/D");
   tree->Branch("diffVertex",&diffVertex,"diffVertex/D");


   //wide jets
   tree->Branch("MJJWide", &MJJWide, "MJJWide/D");
   tree->Branch("DeltaEtaJJWide", &DeltaEtaJJWide, "DeltaEtaJJWide/D");
   tree->Branch("DeltaPhiJJWide", &DeltaPhiJJWide, "DeltaPhiJJWide/D");
   tree->Branch("RJJWide", &RJJWide, "RJJWide/D");
   tree->Branch("wj1Pt", &wj1Pt, "wj1Pt/D");
   tree->Branch("wj2Pt", &wj2Pt, "wj2Pt/D");
   tree->Branch("nGV", &nGV, "nGV/I");
   tree->Branch("nGV_sel", &nGV_sel, "nGV_sel/I");
   
   // histograms

   h_eventSelection = fs->make<TH1F>("eventSelection", " " , 3,0,3);
   h_eventSelection->GetXaxis()->SetBinLabel(1,"All"); 
   h_eventSelection->GetXaxis()->SetBinLabel(2,"jets"); 
   h_eventSelection->GetXaxis()->SetBinLabel(3,"tracks");

   h_pv_multiplicity  = fs->make<TH1F>( "h_pv_multiplicity",  "h_pv_multiplicity", 50,   -0.5,   50.5 );
   h_pv_x             = fs->make<TH1F>( "h_pv_x",             "h_pv_x", 100,   -0.5,  0.5 );
   h_pv_y             = fs->make<TH1F>( "h_pv_y",             "h_pv_y", 100,   -0.5,  0.5 );
   h_pv_z             = fs->make<TH1F>( "h_pv_z",             "h_pv_z", 120,   -30.,  30. );

}


AnalyzerPPS::~AnalyzerPPS()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


/////////////////////////////////////////
//     Variables Initialization        //
/////////////////////////////////////////
void AnalyzerPPS::Init()
{

// wrt MCGenInfo
  Pthat = -1.;
  GeneratorWeight = 1.;
  runNumber = -1;
  eventNumber = -1;
//wrt GenCollection
  GenProtonVectorInfo.clear();
  protonLorentzVector.clear(); 
  Genjetsize = -999;
  GenMxx = -999.;
  GenXi0 = -999.;
  GenXi1 = -999.;
  GenMx = -999.;
  GenMjj = -999.;
  GenRjj = -999.;
  GenJetsVector.clear();
  GenJetsVector_pt.clear();
  GenJetsVector_eta.clear();
  GenJetsVector_phi.clear();
// wrt PPSGenInfo
  PPS_Gen_xiARMF = -999.;
  PPS_Gen_tARMF  = -999.; 
  PPS_Gen_etaARMF= -999.; 
  PPS_Gen_phiARMF= -999.;
  PPS_Gen_xiARMB = -999.; 
  PPS_Gen_tARMB  = -999.; 
  PPS_Gen_etaARMB= -999.;
  PPS_Gen_phiARMB= -999.;
  PPSGenVertexVectorX.clear();
  PPSGenVertexVectorY.clear();
  PPSGenVertexVectorZ.clear();
// wrt FIllCollection
  VertexVector.clear();
  nVertex = 0;
  JetsVector.clear();
  jetsize = -999;
  PPSRecoVertexSize = -999;
  PPSRecoVertexVectorX.clear();
  PPSRecoVertexVectorY.clear();
  PPSRecoVertexVectorZ.clear();
  xiPPSARMF.clear();
  tPPSARMF.clear();
  etaPPSARMF.clear();
  phiPPSARMF.clear();
  thetaXPPSARMF.clear();
  thetaYPPSARMF.clear();
  tofPPSARMF.clear();
  xDet1PPSARMF.clear();
  yDet1PPSARMF.clear();
  xDet2PPSARMF.clear();
  yDet2PPSARMF.clear();
  xTofPPSARMF.clear();
  yTofPPSARMF .clear();
  xiPPSARMB.clear();
  tPPSARMB.clear();
  etaPPSARMB.clear();
  phiPPSARMB.clear();
  thetaXPPSARMB.clear();
  thetaYPPSARMB.clear();
  tofPPSARMB.clear();
  xDet1PPSARMB.clear();
  yDet1PPSARMB.clear();
  xDet2PPSARMB.clear();
  yDet2PPSARMB.clear();
  xTofPPSARMB.clear();
  yTofPPSARMB .clear();
  Mx = -999.;
  VertexVectorX.clear();
  VertexVectorY.clear();
  VertexVectorZ.clear();
  PPSCMSVertex.clear();
  diffVertex = -999.;
  MJJWide = -999.;
  DeltaEtaJJWide = -999.;
  DeltaPhiJJWide = -999.;
  nGV =0;
  nGV_sel=0;
  
  // wrt Sortingobjects
  JetsVector_pt.clear();
  JetsVector_eta.clear();
  JetsVector_phi.clear();
  Mjj= -999.;
  Rjj= -999.;
  yRapidity = -999.;
  DeltaEta12 = -999.;
  DeltaPhi12 = -999.;

  // Jets associate vertex
  GoldenVertexZ = -999.;
  MinimumDistance.clear();
  MinDistance = -1.;
  MaxDistance = -1.;

  
} // end of Init function





/////////////////////////////////////////
//             MC Gen Info             //
/////////////////////////////////////////
void AnalyzerPPS::MCGenInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  unsigned int eventNumber_ = iEvent.id().event();
  unsigned int runNumber_ = iEvent.id().run();
  eventNumber = eventNumber_;
  runNumber = runNumber_;
  if( runWithWeightGen_ ){
    edm::Handle<GenEventInfoProduct> genEventInfoH;
    iEvent.getByLabel("generator", genEventInfoH);
    Pthat = genEventInfoH->binningValues()[0] ;
    GeneratorWeight= genEventInfoH->weight() ;
  } else {
    Pthat= -1. ;
    GeneratorWeight= 1. ;
  }
} // end of MCGenInfo function





/////////////////////////////////////////
//         Gen Collection              //
/////////////////////////////////////////
void AnalyzerPPS::GenCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //Accessing Gen Particle Info
   edm::Handle<reco::GenParticleCollection> genParticle;
   iEvent.getByToken(genParticleToken, genParticle);

   math::XYZTLorentzVector allGenParticles(0.,0.,0.,0.);

   int gensize = genParticle->size();
   int itGen;

   if(genParticle->size()>0){
     for(itGen=0; itGen < gensize; ++itGen){
       const reco::GenParticle* genAll = &((*genParticle)[itGen]);  
       if (genAll->pdgId() == 2212) GenProtonVectorInfo.push_back(genAll);
     }
   }
  
   // Saving One or Two Leading Protons
   if (GenProtonVectorInfo.size()==1){
     protonLorentzVector.push_back(GenProtonVectorInfo[0]->p4());
   }
   if (GenProtonVectorInfo.size()>1){
     protonLorentzVector.push_back(GenProtonVectorInfo[0]->p4());
     protonLorentzVector.push_back(GenProtonVectorInfo[1]->p4());
     GenXi0 = (1-GenProtonVectorInfo[0]->pz()/6500);
     GenXi1 = (1+GenProtonVectorInfo[1]->pz()/6500);
     GenMx = sqrt(GenXi0*GenXi1)*13000;
     // calcular o GenMx usando o xi dos protons!!!
     // see https://arxiv.org/pdf/1006.1289.pdf (page 4)
   }

   //Gen Jets Information
   edm::Handle<reco::GenJetCollection> genjets;
   iEvent.getByToken(genjetsToken,genjets);

   //int Genjetsize = genjets->size();
   Genjetsize = genjets->size();
   int itGenJets; 

   if(genjets->size()>0){
     for(itGenJets=0; itGenJets < Genjetsize; ++itGenJets){
       const reco::GenJet* GenjetAll = &((*genjets)[itGenJets]);
       GenJetsVector.push_back(GenjetAll);
     }
   }

   // Ordering Jets by pT and Fill Jet Vectors
   if(GenJetsVector.size()>0){
     const int GenJetsVectorSize = (int) GenJetsVector.size();
     int *sortGenJetsVector= new int[GenJetsVectorSize];
     double *genvjets = new double[GenJetsVectorSize];
     for (int i=0; i< GenJetsVectorSize; i++) {
       genvjets[i] = GenJetsVector[i]->pt();
     }
     TMath::Sort(GenJetsVectorSize, genvjets, sortGenJetsVector, true);
     for (unsigned int i=0;i<GenJetsVector.size();i++){
       GenJetsVector_pt.push_back(GenJetsVector[sortGenJetsVector[i]]->pt());
       GenJetsVector_eta.push_back(GenJetsVector[sortGenJetsVector[i]]->eta());
       GenJetsVector_phi.push_back(GenJetsVector[sortGenJetsVector[i]]->phi());
     }
   }

   // Gen Mjj and Rjj

   if(Genjetsize < 2) return;
   const reco::GenJet* genJet1 = &(*genjets)[0];
   const reco::GenJet* genJet2 = &(*genjets)[1];
   if(genJet1&&genJet2){
     math::XYZTLorentzVector dijetGenSystem(0.,0.,0.,0.);
     dijetGenSystem += genJet1->p4();
     dijetGenSystem += genJet2->p4();
     double massGen = dijetGenSystem.M();
     GenMjj=massGen;
   }
   if ( (GenMjj > 0.) && (GenMx > 0.)) GenRjj = GenMjj/GenMx; 

} // end of GenCollections function




/////////////////////////////////////////
//           GenPPS Info               //
/////////////////////////////////////////
void AnalyzerPPS::GenPPSInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm; 

  edm::Handle<PPSSpectrometer<PPSGen> > ppsGen;
  iEvent.getByToken(ppsGenToken, ppsGen);

  //Arm F 
  for(size_t i=0;i<ppsGen->ArmF.genParticles.size();++i) {
	PPS_Gen_xiARMF  = ppsGen->ArmF.genParticles.at(i).xi;
	PPS_Gen_tARMF   = ppsGen->ArmF.genParticles.at(i).t;
	PPS_Gen_etaARMF = ppsGen->ArmF.genParticles.at(i).eta;
	PPS_Gen_phiARMF = ppsGen->ArmF.genParticles.at(i).phi;
  }

  //Arm B
  for(size_t i=0;i<ppsGen->ArmB.genParticles.size();++i) {
	PPS_Gen_xiARMB  = ppsGen->ArmB.genParticles.at(i).xi;
	PPS_Gen_tARMB   = ppsGen->ArmB.genParticles.at(i).t;
	PPS_Gen_etaARMB = ppsGen->ArmB.genParticles.at(i).eta;
	PPS_Gen_phiARMB = ppsGen->ArmB.genParticles.at(i).phi;
  }

  int ppsGenVtx = ppsGen->Vertices->size();

  if (ppsGenVtx > 0){
    for (int i=0; i<ppsGenVtx; i++){
      PPSGenVertexVectorX.push_back(ppsGen->Vertices->at(i).x);
      PPSGenVertexVectorY.push_back(ppsGen->Vertices->at(i).y);
      PPSGenVertexVectorZ.push_back(ppsGen->Vertices->at(i).z);
    }
  }


} // enf of GenPPSInfo


/////////////////////////////////////////
//         Fill Collection             //
/////////////////////////////////////////
void AnalyzerPPS::FillCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  h_eventSelection->Fill("All",GeneratorWeight); 

  // Fill CMS Vertex
  Handle<edm::View<reco::Vertex> > vertex;
  iEvent.getByToken(vertexToken, vertex);

  int vertexsize = vertex->size();
  int itVertex;

  if(vertex.isValid()){
    for(itVertex=0; itVertex < vertexsize; ++itVertex){
      const reco::Vertex* vertexAll = &((*vertex)[itVertex]);
      if (vertexAll->isValid()==0) continue; 
      VertexVector.push_back(vertexAll);
      VertexVectorX.push_back(VertexVector[itVertex]->x());
      VertexVectorY.push_back(VertexVector[itVertex]->y());
      VertexVectorZ.push_back(VertexVector[itVertex]->z());
    }
  }

  nVertex = VertexVector.size();

  // Fill Jets no miniAOD ja estao ordenados em Pt decrescente
  Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetsToken, jets);

  jetsize = jets->size();
  int itJets;

  if(jets->size()>0){
    for(itJets=0; itJets < jetsize; ++itJets){
      const pat::Jet* jetAll = &((*jets)[itJets]);
      JetsVector.push_back(jetAll);
      //cout << "Jet [" << itJets << "]  =  "<< jets->at(itJets).pt() << "  GeV" << endl;
    }
  }

  // Fill WideJets
  edm::Handle<std::vector<math::PtEtaPhiMLorentzVector>> wjets;
  iEvent.getByToken(widejetsToken, wjets);
 
  if( wjets->size() == 2 ) {

    wj1.SetPtEtaPhiM(
      wjets->at(0).pt() , wjets->at(0).eta(),
      wjets->at(0).phi(), wjets->at(0).mass()
    );

    wj2.SetPtEtaPhiM(
      wjets->at(1).pt() , wjets->at(1).eta(),
      wjets->at(1).phi(), wjets->at(1).mass()
    );
    
    wdijet = wj1 + wj2;
    MJJWide = wdijet.M();
    DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
    DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));
    wj1Pt = wj1.Pt();
    wj2Pt = wj2.Pt();
    
    cout << "Mjj = "<<MJJWide<<"   *** DeltaEtaJJWide = "<< DeltaEtaJJWide << "   *** DeltaPhiJJWide = "<<DeltaPhiJJWide <<endl;
    
  }
  
  
  // Fill PPS Info
  edm::Handle<PPSSpectrometer<PPSReco> > ppsReco;
  iEvent.getByToken(ppsRecoToken, ppsReco);

  // CT-PPS reconstructed vertex
  PPSRecoVertexSize = ppsReco->Vertices->size();
  
/* decisão quando houver mais de um Golden Vertice
  // aux var
  double Mx_temp=0, RJJWide_temp=0;
  
  for (size_t i=0; i<ppsReco->Vertices->size(); ++i){
	if(ppsReco->Vertices->at(i).Flag){
	    nGV++;
	    if (nGV>1){
	      Mx_temp=EBeam_*TMath::Sqrt((ppsReco->ArmF.Tracks.at(ppsReco->Vertices->at(i).idxTrkF).xi)*(ppsReco->ArmB.Tracks.at(ppsReco->Vertices->at(i).idxTrkB).xi));
	      RJJWide_temp = MJJWide/Mx;
	      if (RJJWide_temp>0.8){
		  Mx = Mx_temp;
		  RJJWide = RJJWide_temp;
		  nGV_sel++;
	      }
	    }
	}
  }
*/

  
  for (size_t i=0; i<ppsReco->Vertices->size(); ++i){
	if(ppsReco->Vertices->at(i).Flag){ //check if it is a golden vertex
	    // calculando a diferença entre o Z do vertice do PPS e do CMS (vertice mais duro)
	    diffVertex = ppsReco->Vertices->at(i).z - VertexVector[0]->z();
	    
	   // if (abs(diffVertex>=PPSres_)){
	      Mx=EBeam_*TMath::Sqrt((ppsReco->ArmF.Tracks.at(ppsReco->Vertices->at(i).idxTrkF).xi)*(ppsReco->ArmB.Tracks.at(ppsReco->Vertices->at(i).idxTrkB).xi));
	   // }
	}
  }
    
  
  for(size_t k=0; k<ppsReco->Vertices->size(); ++k) {
	PPSRecoVertexVectorX.push_back(ppsReco->Vertices->at(k).x);
	PPSRecoVertexVectorY.push_back(ppsReco->Vertices->at(k).y);
	PPSRecoVertexVectorZ.push_back(ppsReco->Vertices->at(k).z);
  }


  //ArmF : xi, t, eta, phi,  thetaX, thetaY, tof
  for(size_t i=0;i<ppsReco->ArmF.Tracks.size();++i) {
	xiPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).xi);
	tPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).t);
	etaPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).eta);
	phiPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).phi);
	thetaXPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).thetaX);
	thetaYPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).thetaY);
	tofPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).ToF.ToF);
	xDet1PPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).Det1.X);
	yDet1PPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).Det1.Y);
	xDet2PPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).Det2.X);
	yDet2PPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).Det2.Y);
	if(ppsReco->ArmF.Tracks.at(i).ToF.ToF!=0){
           xTofPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).ToF.X);
	   yTofPPSARMF.push_back(ppsReco->ArmF.Tracks.at(i).ToF.Y); 
        }
  }

  //ARMB : xi, t, eta, phi,  thetaX, thetaY, tof
  for(size_t i=0;i<ppsReco->ArmB.Tracks.size();++i) {
	xiPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).xi);
	tPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).t);
	etaPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).eta);
	phiPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).phi);
	thetaXPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).thetaX);
	thetaYPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).thetaY);
	tofPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).ToF.ToF);
	xDet1PPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).Det1.X);
	yDet1PPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).Det1.Y);
	xDet2PPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).Det2.X);
	yDet2PPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).Det2.Y);

	// PPS Mx (diffractive mass see http://www.cbpf.br/~gilvan/Proposal/node8.html)
//	if((PPSRecoVertexSize>0)){   
//	   Mx = EBeam_*TMath::Sqrt(xiPPSARMF[0]*xiPPSARMB[0]);
//	}
	if(ppsReco->ArmB.Tracks.at(i).ToF.ToF!=0){
           xTofPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).ToF.X);
	   yTofPPSARMB.push_back(ppsReco->ArmB.Tracks.at(i).ToF.Y); 
        }
  }


  // Ratio MJJWide/Mx
  if(Mx !=0) RJJWide = MJJWide/Mx;

} // enf of FillCollections




/////////////////////////////////////////
//          SortingObjects             //
/////////////////////////////////////////

void AnalyzerPPS::SortingObjects(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  // Ordering Jets by pT and Fill Jet Vectors
  if(JetsVector.size()>0){
    const int JetsVectorSize = (int) JetsVector.size();
    int *sortJetsVector= new int[JetsVectorSize];
    double *vjets = new double[JetsVectorSize];

    for (int i=0; i<JetsVectorSize; i++) {
      vjets[i] = JetsVector[i]->pt();
    }

    TMath::Sort(JetsVectorSize, vjets, sortJetsVector, true);

    for (unsigned int i=0;i<JetsVector.size();i++){
      JetsVector_pt.push_back(JetsVector[sortJetsVector[i]]->pt());
      JetsVector_eta.push_back(JetsVector[sortJetsVector[i]]->eta());
      JetsVector_phi.push_back(JetsVector[sortJetsVector[i]]->phi());
    }

    if(JetsVector.size()>1){
      math::XYZTLorentzVector dijetSystem(0.,0.,0.,0.);
      dijetSystem += JetsVector[sortJetsVector[0]]->p4();
      dijetSystem += JetsVector[sortJetsVector[1]]->p4();
      Mjj = dijetSystem.M();
      yRapidity = dijetSystem.Rapidity();
    }
    
    // Delta phi e Delta eta entre os jatos 12 e 13
    if (JetsVector.size()==2){
      DeltaPhi12 = (JetsVector[sortJetsVector[0]]->phi())-(JetsVector[sortJetsVector[1]]->phi());
      DeltaEta12 = (JetsVector[sortJetsVector[0]]->eta())-(JetsVector[sortJetsVector[1]]->eta());
    }

    if ( (Mjj > 0.) && (Mx > 0.)){
      Rjj = Mjj/Mx;
    }

  }
  
   
} // end of  SortingObjects





/////////////////////////////////////////
//      AssociateJetsWithVertex        //
/////////////////////////////////////////



void AnalyzerPPS::AssociateJetsWithVertex(const edm::Event& iEvent, const edm::EventSetup& iSetup){

/*
  if (indexGold != -999) { 
     GoldenVertexZ = VertexVector[indexGold]->z();
	for (unsigned int i=0;i<TracksVector.size();i++){
	    MinimumDistance.push_back(fabs(VertexVector[indexGold]->z() - TracksVector[i]->vertex().Z()));
 	}
  }  


  const int minVectorSize = (int) MinimumDistance.size();
  int *sortMinVector= new int[minVectorSize];
  double *vmin = new double[minVectorSize];

  for (int i=0; i<minVectorSize; i++) {
    vmin[i]=MinimumDistance[i];
  }
  TMath::Sort(minVectorSize, vmin, sortMinVector, false);
  if(MinimumDistance.size()>0) {
    //cout << "Minimum Distance | Tracks(z)-Vertex(z) |: " << MinimumDistance[sortMinVector[0]] << " cm" << endl;
    //cout << "Maximum Distance | Tracks(z)-Vertex(z) |: " << MinimumDistance[sortMinVector[MinimumDistance.size()-1]] << " cm\n" << endl;
    MinDistance = MinimumDistance[sortMinVector[0]];
    MaxDistance = MinimumDistance[sortMinVector[MinimumDistance.size()-1]];
  }

*/

/*
  if(JetsVector.size()>0 && VertexVector.size()>0){
    h_eventSelection->Fill("jets",GeneratorWeight); 
    math::XYZVector coordleadingjet(JetVertex[0].X(),JetVertex[0].Y(),JetVertex[0].Z());
    JetsSamePosition.push_back(coordleadingjet);
    JetsSameVector_pt.push_back(JetsVector[0]->pt());
    JetsSameVector_eta.push_back(JetsVector[0]->eta());
    JetsSameVector_phi.push_back(JetsVector[0]->phi());
    JetsSameVector_p4.push_back(JetsVector[0]->p4());
    JetsSameVertex_x = JetVertex[0].X();
    JetsSameVertex_y = JetVertex[0].Y();
    JetsSameVertex_z = JetVertex[0].Z();

    for(unsigned int i=1;i<JetsVector.size();i++){
      DistanceBetweenJets.push_back(fabs(JetVertex[0].Z()-JetVertex[i].Z()));      
      if(fabs(JetVertex[0].Z()-JetVertex[i].Z()) < cmsVertexResolution_){
	math::XYZVector coordjets(JetVertex[i].X(),JetVertex[i].Y(),JetVertex[i].Z());
	JetsSamePosition.push_back(coordjets);
	JetsSameVector_pt.push_back(JetsVector[i]->pt());
	JetsSameVector_eta.push_back(JetsVector[i]->eta());
	JetsSameVector_phi.push_back(JetsVector[i]->phi());
	JetsSameVector_p4.push_back(JetsVector[i]->p4());
      }
    }
  }
*/
} // end of AssociateJetsWithVertex



// ------------ method called for each event  ------------
void AnalyzerPPS::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;


////////   Clean Variables   /////////////
   Init(); 

////////     MC Gen Info     /////////////
//  MCGenInfo(iEvent, iSetup);

////////   Gen Collection    /////////////
//  GenCollections(iEvent, iSetup);

////////    PPS Gen Info     /////////////
//  GenPPSInfo(iEvent, iSetup);

////////   Fill Collectios   /////////////
  FillCollections(iEvent, iSetup);

////////   Sorting Objects   /////////////
  SortingObjects(iEvent, iSetup); 

////// Associate Jets With Vertex ////////
//  AssociateJetsWithVertex(iEvent, iSetup);

////////        Fill Tree      ///////////
//////// End Loop on Each Evts ///////////

   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
AnalyzerPPS::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnalyzerPPS::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzerPPS::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzerPPS);
