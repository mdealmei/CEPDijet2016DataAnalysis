//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 21 18:55:48 2016 by ROOT version 5.34/24
// from TTree PPS/PPS
// found on file: AnalyzerPPS_PU.root
//////////////////////////////////////////////////////////

#ifndef makeHistos_h
#define makeHistos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <Math/GenVector/PxPyPzE4D.h>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxProtonsP4 = 2;

class makeHistos {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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

   makeHistos(TTree *tree=0);
   virtual ~makeHistos();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makeHistos_cxx
makeHistos::makeHistos(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("AnalyzerPPS_PU.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("AnalyzerPPS_PU.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("AnalyzerPPS_PU.root:/demo");
      dir->GetObject("PPS",tree);

   }
   Init(tree);
}

makeHistos::~makeHistos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t makeHistos::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makeHistos::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void makeHistos::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   GenJetsPt = 0;
   GenJetsEta = 0;
   GenJetsPhi = 0;
   PPSGenVertexVector_x = 0;
   PPSGenVertexVector_y = 0;
   PPSGenVertexVector_z = 0;
   PPSRecoVertexVector_x = 0;
   PPSRecoVertexVector_y = 0;
   PPSRecoVertexVector_z = 0;
   xiPPSARMF = 0;
   tPPSARMF = 0;
   etaPPSARMF = 0;
   phiPPSARMF = 0;
   thetaXPPSARMF = 0;
   thetaYPPSARMF = 0;
   tofPPSARMF = 0;
   xDet1PPSARMF = 0;
   yDet1PPSARMF = 0;
   xDet2PPSARMF = 0;
   yDet2PPSARMF = 0;
   xTofPPSARMF = 0;
   yTofPPSARMF = 0;
   xiPPSARMB = 0;
   tPPSARMB = 0;
   etaPPSARMB = 0;
   phiPPSARMB = 0;
   thetaXPPSARMB = 0;
   thetaYPPSARMB = 0;
   tofPPSARMB = 0;
   xDet1PPSARMB = 0;
   yDet1PPSARMB = 0;
   xDet2PPSARMB = 0;
   yDet2PPSARMB = 0;
   xTofPPSARMB = 0;
   yTofPPSARMB = 0;
   VertexVectorX = 0;
   VertexVectorY = 0;
   VertexVectorZ = 0;
   JetsPt = 0;
   JetsEta = 0;
   JetsPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pv_multiplicity", &pv_multiplicity, &b_pv_multiplicity);
   fChain->SetBranchAddress("GeneratorWeight", &GeneratorWeight, &b_GeneratorWeight);
   fChain->SetBranchAddress("ProtonsP4", &ProtonsP4_, &b_ProtonsP4_);
   fChain->SetBranchAddress("ProtonsP4.fCoordinates.fX", ProtonsP4_fCoordinates_fX, &b_ProtonsP4_fCoordinates_fX);
   fChain->SetBranchAddress("ProtonsP4.fCoordinates.fY", ProtonsP4_fCoordinates_fY, &b_ProtonsP4_fCoordinates_fY);
   fChain->SetBranchAddress("ProtonsP4.fCoordinates.fZ", ProtonsP4_fCoordinates_fZ, &b_ProtonsP4_fCoordinates_fZ);
   fChain->SetBranchAddress("ProtonsP4.fCoordinates.fT", ProtonsP4_fCoordinates_fT, &b_ProtonsP4_fCoordinates_fT);
   fChain->SetBranchAddress("GenMxx", &GenMxx, &b_GenMxx);
   fChain->SetBranchAddress("GenJetsMultiplicity", &GenJetsMultiplicity, &b_GenJetsMultiplicity);
   fChain->SetBranchAddress("GenJetsPt", &GenJetsPt, &b_GenJetsPt);
   fChain->SetBranchAddress("GenJetsEta", &GenJetsEta, &b_GenJetsEta);
   fChain->SetBranchAddress("GenJetsPhi", &GenJetsPhi, &b_GenJetsPhi);
   fChain->SetBranchAddress("GenMjj", &GenMjj, &b_GenMjj);
   fChain->SetBranchAddress("GenRjj", &GenRjj, &b_GenRjj);
   fChain->SetBranchAddress("PPS_Gen_xiARMF", &PPS_Gen_xiARMF, &b_PPS_Gen_xiARMF);
   fChain->SetBranchAddress("PPS_Gen_tARMF", &PPS_Gen_tARMF, &b_PPS_Gen_tARMF);
   fChain->SetBranchAddress("PPS_Gen_etaARMF", &PPS_Gen_etaARMF, &b_PPS_Gen_etaARMF);
   fChain->SetBranchAddress("PPS_Gen_phiARMF", &PPS_Gen_phiARMF, &b_PPS_Gen_phiARMF);
   fChain->SetBranchAddress("PPS_Gen_xiARMB", &PPS_Gen_xiARMB, &b_PPS_Gen_xiARMB);
   fChain->SetBranchAddress("PPS_Gen_tARMB", &PPS_Gen_tARMB, &b_PPS_Gen_tARMB);
   fChain->SetBranchAddress("PPS_Gen_etaARMB", &PPS_Gen_etaARMB, &b_PPS_Gen_etaARMB);
   fChain->SetBranchAddress("PPS_Gen_phiARMB", &PPS_Gen_phiARMB, &b_PPS_Gen_phiARMB);
   fChain->SetBranchAddress("PPSGenVertexVector_x", &PPSGenVertexVector_x, &b_PPSGenVertexVector_x);
   fChain->SetBranchAddress("PPSGenVertexVector_y", &PPSGenVertexVector_y, &b_PPSGenVertexVector_y);
   fChain->SetBranchAddress("PPSGenVertexVector_z", &PPSGenVertexVector_z, &b_PPSGenVertexVector_z);
   fChain->SetBranchAddress("nVertex", &nVertex, &b_nVertex);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("Mpf", &Mpf, &b_Mpf);
   fChain->SetBranchAddress("PPSRecoVertexSize", &PPSRecoVertexSize, &b_PPSRecoVertexSize);
   fChain->SetBranchAddress("PPSRecoVertexVector_x", &PPSRecoVertexVector_x, &b_PPSRecoVertexVector_x);
   fChain->SetBranchAddress("PPSRecoVertexVector_y", &PPSRecoVertexVector_y, &b_PPSRecoVertexVector_y);
   fChain->SetBranchAddress("PPSRecoVertexVector_z", &PPSRecoVertexVector_z, &b_PPSRecoVertexVector_z);
   fChain->SetBranchAddress("xiPPSARMF", &xiPPSARMF, &b_xiPPSARMF);
   fChain->SetBranchAddress("tPPSARMF", &tPPSARMF, &b_tPPSARMF);
   fChain->SetBranchAddress("etaPPSARMF", &etaPPSARMF, &b_etaPPSARMF);
   fChain->SetBranchAddress("phiPPSARMF", &phiPPSARMF, &b_phiPPSARMF);
   fChain->SetBranchAddress("thetaXPPSARMF", &thetaXPPSARMF, &b_thetaXPPSARMF);
   fChain->SetBranchAddress("thetaYPPSARMF", &thetaYPPSARMF, &b_thetaYPPSARMF);
   fChain->SetBranchAddress("tofPPSARMF", &tofPPSARMF, &b_tofPPSARMF);
   fChain->SetBranchAddress("xDet1PPSARMF", &xDet1PPSARMF, &b_xDet1PPSARMF);
   fChain->SetBranchAddress("yDet1PPSARMF", &yDet1PPSARMF, &b_yDet1PPSARMF);
   fChain->SetBranchAddress("xDet2PPSARMF", &xDet2PPSARMF, &b_xDet2PPSARMF);
   fChain->SetBranchAddress("yDet2PPSARMF", &yDet2PPSARMF, &b_yDet2PPSARMF);
   fChain->SetBranchAddress("xTofPPSARMF", &xTofPPSARMF, &b_xTofPPSARMF);
   fChain->SetBranchAddress("yTofPPSARMF", &yTofPPSARMF, &b_yTofPPSARMF);
   fChain->SetBranchAddress("xiPPSARMB", &xiPPSARMB, &b_xiPPSARMB);
   fChain->SetBranchAddress("tPPSARMB", &tPPSARMB, &b_tPPSARMB);
   fChain->SetBranchAddress("etaPPSARMB", &etaPPSARMB, &b_etaPPSARMB);
   fChain->SetBranchAddress("phiPPSARMB", &phiPPSARMB, &b_phiPPSARMB);
   fChain->SetBranchAddress("thetaXPPSARMB", &thetaXPPSARMB, &b_thetaXPPSARMB);
   fChain->SetBranchAddress("thetaYPPSARMB", &thetaYPPSARMB, &b_thetaYPPSARMB);
   fChain->SetBranchAddress("tofPPSARMB", &tofPPSARMB, &b_tofPPSARMB);
   fChain->SetBranchAddress("xDet1PPSARMB", &xDet1PPSARMB, &b_xDet1PPSARMB);
   fChain->SetBranchAddress("yDet1PPSARMB", &yDet1PPSARMB, &b_yDet1PPSARMB);
   fChain->SetBranchAddress("xDet2PPSARMB", &xDet2PPSARMB, &b_xDet2PPSARMB);
   fChain->SetBranchAddress("yDet2PPSARMB", &yDet2PPSARMB, &b_yDet2PPSARMB);
   fChain->SetBranchAddress("xTofPPSARMB", &xTofPPSARMB, &b_xTofPPSARMB);
   fChain->SetBranchAddress("yTofPPSARMB", &yTofPPSARMB, &b_yTofPPSARMB);
   fChain->SetBranchAddress("Mx", &Mx, &b_Mx);
   fChain->SetBranchAddress("VertexVectorX", &VertexVectorX, &b_VertexVectorX);
   fChain->SetBranchAddress("VertexVectorY", &VertexVectorY, &b_VertexVectorY);
   fChain->SetBranchAddress("VertexVectorZ", &VertexVectorZ, &b_VertexVectorZ);
   fChain->SetBranchAddress("JetsPt", &JetsPt, &b_JetsPt);
   fChain->SetBranchAddress("JetsEta", &JetsEta, &b_JetsEta);
   fChain->SetBranchAddress("JetsPhi", &JetsPhi, &b_JetsPhi);
   fChain->SetBranchAddress("Mjj", &Mjj, &b_Mjj);
   fChain->SetBranchAddress("Rjj", &Rjj, &b_Rjj);
   fChain->SetBranchAddress("yRapidity", &yRapidity, &b_yRapidity);
   fChain->SetBranchAddress("GoldenVertexZ", &GoldenVertexZ, &b_GoldenVertexZ);
   fChain->SetBranchAddress("MinDistance", &MinDistance, &b_MinDistance);
   fChain->SetBranchAddress("MaxDistance", &MaxDistance, &b_MaxDistance);
}

void makeHistos::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

#endif // #ifdef makeHistos_cxx
