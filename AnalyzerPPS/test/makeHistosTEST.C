#include <stdio.h>      /* printf */
#include <math.h>       /* sqrt */
//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>
#include <TVectorT.h>
//
////STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
//#include <vector>
#include <map>
#include <cmath>
#define PI 3.141592653589793
using namespace std;

void makeHistosTEST(){

  //------------ FWLite libraries ------
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  gROOT->ProcessLine("#include<vector>");
  gROOT->ProcessLine(".exception"); 
       
  // input and output files
  TFile *in  = TFile::Open("AnalyzerPPS_PU.root");
  TFile *out = new TFile("makeHistos.root", "RECREATE");

  // accessing tree
  TTree *tr = (TTree*)in->Get("demo/PPS");

   // Declaration of leaf types
   Double_t        pv_multiplicity;
   Double_t        GeneratorWeight;
   Int_t           ProtonsP4_;
   Double_t        ProtonsP4_fCoordinates_fX[kMaxProtonsP4];   //[ProtonsP4_]
   Double_t        ProtonsP4_fCoordinates_fY[kMaxProtonsP4];   //[ProtonsP4_]
   Double_t        ProtonsP4_fCoordinates_fZ[kMaxProtonsP4];   //[ProtonsP4_]
   Double_t        ProtonsP4_fCoordinates_fT[kMaxProtonsP4];   //[ProtonsP4_]
   Double_t        GenMxx;
   Int_t           GenJetsMultiplicity;
   vector<double>  *GenJetsPt;
   vector<double>  *GenJetsEta;
   vector<double>  *GenJetsPhi;
   Double_t        GenMjj;
   Double_t        GenRjj;
   Double_t        PPS_Gen_xiARMF;
   Double_t        PPS_Gen_tARMF;
   Double_t        PPS_Gen_etaARMF;
   Double_t        PPS_Gen_phiARMF;
   Double_t        PPS_Gen_xiARMB;
   Double_t        PPS_Gen_tARMB;
   Double_t        PPS_Gen_etaARMB;
   Double_t        PPS_Gen_phiARMB;
   vector<double>  *PPSGenVertexVector_x;
   vector<double>  *PPSGenVertexVector_y;
   vector<double>  *PPSGenVertexVector_z;
   Int_t           nVertex;
   Int_t           nTracks;
   Double_t        Mpf;
   Int_t           PPSRecoVertexSize;
   vector<double>  *PPSRecoVertexVector_x;
   vector<double>  *PPSRecoVertexVector_y;
   vector<double>  *PPSRecoVertexVector_z;
   vector<double>  *xiPPSARMF;
   vector<double>  *tPPSARMF;
   vector<double>  *etaPPSARMF;
   vector<double>  *phiPPSARMF;
   vector<double>  *thetaXPPSARMF;
   vector<double>  *thetaYPPSARMF;
   vector<double>  *tofPPSARMF;
   vector<double>  *xDet1PPSARMF;
   vector<double>  *yDet1PPSARMF;
   vector<double>  *xDet2PPSARMF;
   vector<double>  *yDet2PPSARMF;
   vector<double>  *xTofPPSARMF;
   vector<double>  *yTofPPSARMF;
   vector<double>  *xiPPSARMB;
   vector<double>  *tPPSARMB;
   vector<double>  *etaPPSARMB;
   vector<double>  *phiPPSARMB;
   vector<double>  *thetaXPPSARMB;
   vector<double>  *thetaYPPSARMB;
   vector<double>  *tofPPSARMB;
   vector<double>  *xDet1PPSARMB;
   vector<double>  *yDet1PPSARMB;
   vector<double>  *xDet2PPSARMB;
   vector<double>  *yDet2PPSARMB;
   vector<double>  *xTofPPSARMB;
   vector<double>  *yTofPPSARMB;
   Double_t        Mx;
   vector<double>  *VertexVectorX;
   vector<double>  *VertexVectorY;
   vector<double>  *VertexVectorZ;
   vector<double>  *JetsPt;
   vector<double>  *JetsEta;
   vector<double>  *JetsPhi;
   Double_t        Mjj;
   Double_t        Rjj;
   Double_t        yRapidity;
   Double_t        GoldenVertexZ;
   Double_t        MinDistance;
   Double_t        MaxDistance;

   // List of branches
   TBranch        *b_pv_multiplicity;   //!
   TBranch        *b_GeneratorWeight;   //!
   TBranch        *b_ProtonsP4_;   //!
   TBranch        *b_ProtonsP4_fCoordinates_fX;   //!
   TBranch        *b_ProtonsP4_fCoordinates_fY;   //!
   TBranch        *b_ProtonsP4_fCoordinates_fZ;   //!
   TBranch        *b_ProtonsP4_fCoordinates_fT;   //!
   TBranch        *b_GenMxx;   //!
   TBranch        *b_GenJetsMultiplicity;   //!
   TBranch        *b_GenJetsPt;   //!
   TBranch        *b_GenJetsEta;   //!
   TBranch        *b_GenJetsPhi;   //!
   TBranch        *b_GenMjj;   //!
   TBranch        *b_GenRjj;   //!
   TBranch        *b_PPS_Gen_xiARMF;   //!
   TBranch        *b_PPS_Gen_tARMF;   //!
   TBranch        *b_PPS_Gen_etaARMF;   //!
   TBranch        *b_PPS_Gen_phiARMF;   //!
   TBranch        *b_PPS_Gen_xiARMB;   //!
   TBranch        *b_PPS_Gen_tARMB;   //!
   TBranch        *b_PPS_Gen_etaARMB;   //!
   TBranch        *b_PPS_Gen_phiARMB;   //!
   TBranch        *b_PPSGenVertexVector_x;   //!
   TBranch        *b_PPSGenVertexVector_y;   //!
   TBranch        *b_PPSGenVertexVector_z;   //!
   TBranch        *b_nVertex;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_Mpf;   //!
   TBranch        *b_PPSRecoVertexSize;   //!
   TBranch        *b_PPSRecoVertexVector_x;   //!
   TBranch        *b_PPSRecoVertexVector_y;   //!
   TBranch        *b_PPSRecoVertexVector_z;   //!
   TBranch        *b_xiPPSARMF;   //!
   TBranch        *b_tPPSARMF;   //!
   TBranch        *b_etaPPSARMF;   //!
   TBranch        *b_phiPPSARMF;   //!
   TBranch        *b_thetaXPPSARMF;   //!
   TBranch        *b_thetaYPPSARMF;   //!
   TBranch        *b_tofPPSARMF;   //!
   TBranch        *b_xDet1PPSARMF;   //!
   TBranch        *b_yDet1PPSARMF;   //!
   TBranch        *b_xDet2PPSARMF;   //!
   TBranch        *b_yDet2PPSARMF;   //!
   TBranch        *b_xTofPPSARMF;   //!
   TBranch        *b_yTofPPSARMF;   //!
   TBranch        *b_xiPPSARMB;   //!
   TBranch        *b_tPPSARMB;   //!
   TBranch        *b_etaPPSARMB;   //!
   TBranch        *b_phiPPSARMB;   //!
   TBranch        *b_thetaXPPSARMB;   //!
   TBranch        *b_thetaYPPSARMB;   //!
   TBranch        *b_tofPPSARMB;   //!
   TBranch        *b_xDet1PPSARMB;   //!
   TBranch        *b_yDet1PPSARMB;   //!
   TBranch        *b_xDet2PPSARMB;   //!
   TBranch        *b_yDet2PPSARMB;   //!
   TBranch        *b_xTofPPSARMB;   //!
   TBranch        *b_yTofPPSARMB;   //!
   TBranch        *b_Mx;   //!
   TBranch        *b_VertexVectorX;   //!
   TBranch        *b_VertexVectorY;   //!
   TBranch        *b_VertexVectorZ;   //!
   TBranch        *b_JetsPt;   //!
   TBranch        *b_JetsEta;   //!
   TBranch        *b_JetsPhi;   //!
   TBranch        *b_Mjj;   //!
   TBranch        *b_Rjj;   //!
   TBranch        *b_yRapidity;   //!
   TBranch        *b_GoldenVertexZ;   //!
   TBranch        *b_MinDistance;   //!
   TBranch        *b_MaxDistance;   //!

   tr->SetBranchAddress("pv_multiplicity", &pv_multiplicity, &b_pv_multiplicity);
   tr->SetBranchAddress("GeneratorWeight", &GeneratorWeight, &b_GeneratorWeight);
   tr->SetBranchAddress("ProtonsP4", &ProtonsP4_, &b_ProtonsP4_);
   tr->SetBranchAddress("ProtonsP4.fCoordinates.fX", ProtonsP4_fCoordinates_fX, &b_ProtonsP4_fCoordinates_fX);
   tr->SetBranchAddress("ProtonsP4.fCoordinates.fY", ProtonsP4_fCoordinates_fY, &b_ProtonsP4_fCoordinates_fY);
   tr->SetBranchAddress("ProtonsP4.fCoordinates.fZ", ProtonsP4_fCoordinates_fZ, &b_ProtonsP4_fCoordinates_fZ);
   tr->SetBranchAddress("ProtonsP4.fCoordinates.fT", ProtonsP4_fCoordinates_fT, &b_ProtonsP4_fCoordinates_fT);
   tr->SetBranchAddress("GenMxx", &GenMxx, &b_GenMxx);
   tr->SetBranchAddress("GenJetsMultiplicity", &GenJetsMultiplicity, &b_GenJetsMultiplicity);
   tr->SetBranchAddress("GenJetsPt", &GenJetsPt, &b_GenJetsPt);
   tr->SetBranchAddress("GenJetsEta", &GenJetsEta, &b_GenJetsEta);
   tr->SetBranchAddress("GenJetsPhi", &GenJetsPhi, &b_GenJetsPhi);
   tr->SetBranchAddress("GenMjj", &GenMjj, &b_GenMjj);
   tr->SetBranchAddress("GenRjj", &GenRjj, &b_GenRjj);
   tr->SetBranchAddress("PPS_Gen_xiARMF", &PPS_Gen_xiARMF, &b_PPS_Gen_xiARMF);
   tr->SetBranchAddress("PPS_Gen_tARMF", &PPS_Gen_tARMF, &b_PPS_Gen_tARMF);
   tr->SetBranchAddress("PPS_Gen_etaARMF", &PPS_Gen_etaARMF, &b_PPS_Gen_etaARMF);
   tr->SetBranchAddress("PPS_Gen_phiARMF", &PPS_Gen_phiARMF, &b_PPS_Gen_phiARMF);
   tr->SetBranchAddress("PPS_Gen_xiARMB", &PPS_Gen_xiARMB, &b_PPS_Gen_xiARMB);
   tr->SetBranchAddress("PPS_Gen_tARMB", &PPS_Gen_tARMB, &b_PPS_Gen_tARMB);
   tr->SetBranchAddress("PPS_Gen_etaARMB", &PPS_Gen_etaARMB, &b_PPS_Gen_etaARMB);
   tr->SetBranchAddress("PPS_Gen_phiARMB", &PPS_Gen_phiARMB, &b_PPS_Gen_phiARMB);
   tr->SetBranchAddress("PPSGenVertexVector_x", &PPSGenVertexVector_x, &b_PPSGenVertexVector_x);
   tr->SetBranchAddress("PPSGenVertexVector_y", &PPSGenVertexVector_y, &b_PPSGenVertexVector_y);
   tr->SetBranchAddress("PPSGenVertexVector_z", &PPSGenVertexVector_z, &b_PPSGenVertexVector_z);
   tr->SetBranchAddress("nVertex", &nVertex, &b_nVertex);
   tr->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   tr->SetBranchAddress("Mpf", &Mpf, &b_Mpf);
   tr->SetBranchAddress("PPSRecoVertexSize", &PPSRecoVertexSize, &b_PPSRecoVertexSize);
   tr->SetBranchAddress("PPSRecoVertexVector_x", &PPSRecoVertexVector_x, &b_PPSRecoVertexVector_x);
   tr->SetBranchAddress("PPSRecoVertexVector_y", &PPSRecoVertexVector_y, &b_PPSRecoVertexVector_y);
   tr->SetBranchAddress("PPSRecoVertexVector_z", &PPSRecoVertexVector_z, &b_PPSRecoVertexVector_z);
   tr->SetBranchAddress("xiPPSARMF", &xiPPSARMF, &b_xiPPSARMF);
   tr->SetBranchAddress("tPPSARMF", &tPPSARMF, &b_tPPSARMF);
   tr->SetBranchAddress("etaPPSARMF", &etaPPSARMF, &b_etaPPSARMF);
   tr->SetBranchAddress("phiPPSARMF", &phiPPSARMF, &b_phiPPSARMF);
   tr->SetBranchAddress("thetaXPPSARMF", &thetaXPPSARMF, &b_thetaXPPSARMF);
   tr->SetBranchAddress("thetaYPPSARMF", &thetaYPPSARMF, &b_thetaYPPSARMF);
   tr->SetBranchAddress("tofPPSARMF", &tofPPSARMF, &b_tofPPSARMF);
   tr->SetBranchAddress("xDet1PPSARMF", &xDet1PPSARMF, &b_xDet1PPSARMF);
   tr->SetBranchAddress("yDet1PPSARMF", &yDet1PPSARMF, &b_yDet1PPSARMF);
   tr->SetBranchAddress("xDet2PPSARMF", &xDet2PPSARMF, &b_xDet2PPSARMF);
   tr->SetBranchAddress("yDet2PPSARMF", &yDet2PPSARMF, &b_yDet2PPSARMF);
   tr->SetBranchAddress("xTofPPSARMF", &xTofPPSARMF, &b_xTofPPSARMF);
   tr->SetBranchAddress("yTofPPSARMF", &yTofPPSARMF, &b_yTofPPSARMF);
   tr->SetBranchAddress("xiPPSARMB", &xiPPSARMB, &b_xiPPSARMB);
   tr->SetBranchAddress("tPPSARMB", &tPPSARMB, &b_tPPSARMB);
   tr->SetBranchAddress("etaPPSARMB", &etaPPSARMB, &b_etaPPSARMB);
   tr->SetBranchAddress("phiPPSARMB", &phiPPSARMB, &b_phiPPSARMB);
   tr->SetBranchAddress("thetaXPPSARMB", &thetaXPPSARMB, &b_thetaXPPSARMB);
   tr->SetBranchAddress("thetaYPPSARMB", &thetaYPPSARMB, &b_thetaYPPSARMB);
   tr->SetBranchAddress("tofPPSARMB", &tofPPSARMB, &b_tofPPSARMB);
   tr->SetBranchAddress("xDet1PPSARMB", &xDet1PPSARMB, &b_xDet1PPSARMB);
   tr->SetBranchAddress("yDet1PPSARMB", &yDet1PPSARMB, &b_yDet1PPSARMB);
   tr->SetBranchAddress("xDet2PPSARMB", &xDet2PPSARMB, &b_xDet2PPSARMB);
   tr->SetBranchAddress("yDet2PPSARMB", &yDet2PPSARMB, &b_yDet2PPSARMB);
   tr->SetBranchAddress("xTofPPSARMB", &xTofPPSARMB, &b_xTofPPSARMB);
   tr->SetBranchAddress("yTofPPSARMB", &yTofPPSARMB, &b_yTofPPSARMB);
   tr->SetBranchAddress("Mx", &Mx, &b_Mx);
   tr->SetBranchAddress("VertexVectorX", &VertexVectorX, &b_VertexVectorX);
   tr->SetBranchAddress("VertexVectorY", &VertexVectorY, &b_VertexVectorY);
   tr->SetBranchAddress("VertexVectorZ", &VertexVectorZ, &b_VertexVectorZ);
   tr->SetBranchAddress("JetsPt", &JetsPt, &b_JetsPt);
   tr->SetBranchAddress("JetsEta", &JetsEta, &b_JetsEta);
   tr->SetBranchAddress("JetsPhi", &JetsPhi, &b_JetsPhi);
   tr->SetBranchAddress("Mjj", &Mjj, &b_Mjj);
   tr->SetBranchAddress("Rjj", &Rjj, &b_Rjj);
   tr->SetBranchAddress("yRapidity", &yRapidity, &b_yRapidity);
   tr->SetBranchAddress("GoldenVertexZ", &GoldenVertexZ, &b_GoldenVertexZ);
   tr->SetBranchAddress("MinDistance", &MinDistance, &b_MinDistance);
   tr->SetBranchAddress("MaxDistance", &MaxDistance, &b_MaxDistance);

    //----------- settings ---------------
    float cross_section =1700.0 ; //[fb]
    float luminosity=100.0;//[fb]^-1
    double lumiweight = (luminosity*cross_section)/NEVENTS;
    cout<<"lumiweight= "<<lumiweight<<endl;
    double mclumiweight = 1.0;

    TH1F *hNJets = new TH1F("NJets","N Jets;  N Jets; N events",100,0,100);
    TH1F *hVertexZCMS = new TH1F("VertexZCMS"," Vertex Z CMS[cm];  Vertex Z [cm]; N events",25,-25.0,25.0);
    TH1F *hVertexZPPS = new TH1F("VertexZPPS","Vertex Z PPS [cm]; Vertex  Z [cm]; N events",25,-25.0,25.0);
    TH2F *hVertexZCMSPPS = new TH2F("VertexZCMSPPS","Vertex Z CMS vs Vertex Z  PPS; Vertex Z CMS [cm]; Vertex Z PPS [cm]",25,-25.0,25.0,25, -25.0,25.0);
    TH1F *hLeadingJetPt = new TH1F("LeadingJetPt","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
    TH1F *hSecondJetPt = new TH1F("SecondJetPt","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
    TH1F *hLeadingJetEta = new TH1F("LeadingJetEta","Leading Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
    TH1F *hSecondJetEta = new TH1F("SecondJetEta","Second Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
    TH1F *hLeadingJetPhi = new TH1F("LeadingJetPhi","Leading Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);
    TH1F *hSecondJetPhi = new TH1F("SecondJetPhi","Second Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);  
    TH1F *hDeltaEtaJets = new TH1F("DeltaEtaJets","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; N events",20,0.0,5.2);
    TH1F *hDeltaPhiJets = new TH1F("DeltaPhiJets","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; N events",20,0.0,1.2*PI);
    TH1F *hMjj = new TH1F( "Mjj" , "Mass_{JJ} Distribution; M_{jj}  [GeV]; N events" , 100, 0., 2000. );
    TH1F *hMx = new TH1F("Mx" , "Mass_{X} Distribution; M_{x}  [GeV]; N events" , 100, 0., 2000. );

    //----------- tree reading -------------------
    unsigned NEntries = tr->GetEntries();
    cout<<"Reading TREE: "<<NEntries<<" events"<<endl;

    int decade = 0;

    for(unsigned i=0;i<NEntries;i++) {

        //----------- progress report -------------
        double progress = 10.0*i/(1.0*NEVENTS);
        Int_t k = TMath::FloorNint(progress);
        if (k > decade) cout<<10*k<<" %"<<endl;
        decade = k;  

        //----------- read the event --------------
        tr->GetEntry(i);
        Njets = JetsPt->size();

        if(Njets<2)continue;



    }


} //end of everything
