// -*- C++ -*-
//
// Package:    CEPDiJetsAnalyzer/fatjets
// Class:      fatjets
// 
/**\class fatjets fatjets.cc CEPDiJetsAnalyzer/fatjets/plugins/fatjets.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Andre De Almeida Pacheco
//         Created:  Sun, 02 Feb 2017 20:16:47 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// Central Detector (CMS)
// Jets
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "TVector3.h"
#include "TLorentzVector.h"

//
// class declaration
//

class fatjets : public edm::stream::EDFilter<> {
   public:
      explicit fatjets(const edm::ParameterSet&);
      ~fatjets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      
      // jets
      edm::EDGetTokenT<pat::JetCollection> jetsToken;
      double wideJetDeltaR_; // Radius parameter for wide jets


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
fatjets::fatjets(const edm::ParameterSet& iConfig):
      jetsToken (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
      wideJetDeltaR_   (iConfig.getParameter<double>("wideJetDeltaR"))
{

  //register your products
  produces<std::vector<math::PtEtaPhiMLorentzVector> >("widejets");

  //now do what ever initialization is needed

}


fatjets::~fatjets()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


/*
Using wide jets: All other jets
with pT > 10 GeV and |η| < 2.5
are added to the closest leading
jet if they are within ΔR<1.1
*/


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
fatjets::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  std::unique_ptr<std::vector<math::PtEtaPhiMLorentzVector> > widejets(new std::vector<math::PtEtaPhiMLorentzVector>); 
  
  // Handle Jets
  Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetsToken, jets);

//  int jetsize = jets->size();
//  cout << "jet multiplicity = " << jetsize << endl;
   
  if(jets->size()<2) return false;

  bool tightJetID=false;
  
  TLorentzVector wj1_tmp;
  TLorentzVector wj2_tmp;
  TLorentzVector wj1;
  TLorentzVector wj2;
  TLorentzVector wdijet;
    
  // Take two leading jets
  TLorentzVector jet1, jet2;

  pat::JetCollection::const_iterator j1 = jets->begin();
  pat::JetCollection::const_iterator j2 = j1; ++j2;     
 
  jet1.SetPtEtaPhiM(j1->pt(),j1->eta(),j1->phi(),j1->mass());
  jet2.SetPtEtaPhiM(j2->pt(),j2->eta(),j2->phi(),j2->mass());  
  
  if (jet1.Pt()<100 || jet2.Pt()<100) return false;  // possivelmente mudar para corte em 100 GeV
  if (abs(jet1.Eta())>2.4 || abs(jet2.Eta())>2.4) return false;
  
//  cout << "j1: Pt = " << jet1.Pt() << "; Eta = " << jet1.Eta() << "; Phi =" << jet1.Phi() << endl;
//  cout << "j2: Pt = " << jet2.Pt() << "; Eta = " << jet2.Eta() << "; Phi =" << jet2.Phi() << endl;
    
  // Create wide jets (radiation recovery algorithm)  
  for(pat::JetCollection::const_iterator it = jets->begin(); it != jets->end(); it++){
  
    // Apply JetID criteria
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_2016
    double NHF = it->neutralHadronEnergyFraction();
    double NEMF = it->neutralEmEnergyFraction();
    double CHF = it->chargedHadronEnergyFraction();
    double CEMF = it->chargedEmEnergyFraction();
    int NumConst = it->chargedMultiplicity()+it->neutralMultiplicity();
    int CHM = it->chargedMultiplicity();
    double eta = it->eta();

    // For |eta|<=2.7 Apply
    if ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7)
      tightJetID = true;
    else
      tightJetID = false;    
   
    //verificar se os 2 leading jets obedecem o JetID antes de passar para os demais
    if((it==j1 ||it==j2) && tightJetID==false) {
//      cout << std::distance(jets->begin(), it) <<"o leading jet not passed on JetID criteria" << endl;
      return false;
    } 
    /*
    else cout << std::distance(jets->begin(), it) <<"o leading jet passed on JetID criteria" << endl;
    
    cout <<"eta  =  " << eta << endl;
    cout <<"NHF  =  " << NHF << endl;
    cout <<"NEMF =  " << NEMF<< endl;
    cout <<"CHF  =  " << CHF << endl;
    cout <<"CEMF =  " << CEMF<< endl;
    cout <<"NumConst =  "<< NumConst << endl;
    cout <<"CHM  =  " << CHM << endl;
    */
    
    TLorentzVector currentJet;
    currentJet.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),it->mass());

    if (currentJet.Pt()<10||abs(currentJet.Eta())>2.4||tightJetID==false) continue;
    // jetID criteria for all jets 
    // https://indico.cern.ch/event/470463/contributions/2153645/attachments/1271951/1886044/Tziaferi_CMS_dijets_13TeV.pdf
    
    double DeltaR1 = currentJet.DeltaR(jet1);
    double DeltaR2 = currentJet.DeltaR(jet2);

    if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_) {
      wj1_tmp += currentJet;    
    } else if(DeltaR2 < wideJetDeltaR_) {
      wj2_tmp += currentJet;
    }
      
      
    //cout << "Jet [" << itJets << "]  Pt =  "<< jets->at(itJets).pt() << "  GeV" << " ***  "<<"eta =  " << jets->at(itJets).eta()<<endl;
    
  }
    
  // Re-order the wide jets in pT 
  if( wj1_tmp.Pt() > wj2_tmp.Pt() ) {
    wj1 = wj1_tmp;
    wj2 = wj2_tmp;     
  }
  else {
    wj1 = wj2_tmp;
    wj2 = wj1_tmp;  
  }

  // Create dijet system
  wdijet = wj1 + wj2;
  
//  cout << "Wide DiJets Mass =  "<< wdijet.M() << " GeV" << endl; 
    
  // Put widejets in the container 
  math::PtEtaPhiMLorentzVector wj1math(wj1.Pt(), wj1.Eta(), wj1.Phi(), wj1.M());
  math::PtEtaPhiMLorentzVector wj2math(wj2.Pt(), wj2.Eta(), wj2.Phi(), wj2.M());
  widejets->push_back( wj1math );
  widejets->push_back( wj2math );
  
  // ## Put objects in the Event
  iEvent.put(std::move(widejets), "widejets");
 
  return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
fatjets::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
fatjets::endStream() {
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
fatjets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(fatjets);