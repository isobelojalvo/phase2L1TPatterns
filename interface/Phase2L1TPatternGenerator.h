/* 
   L1Trigger/phase2L1TPatterns/interface/Phase2L1TPatternGenerator.h

 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

//track trigger data formats
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"


#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <TLorentzVector.h>
#include <memory>
#include <math.h>
#include <vector>
#include <list>
#include <algorithm>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "L1Trigger/phase2L1TPatterns/interface/triggerGeometryTools.hh"

#include <fstream>

using namespace edm;
using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
//
// class declaration
//

class Phase2L1TPatternGenerator : public edm::EDAnalyzer {
   public:
      explicit Phase2L1TPatternGenerator(const edm::ParameterSet&);
      ~Phase2L1TPatternGenerator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   protected:
      void printNZeros(ofstream &file, int nZeros);
      void printEGTau(ofstream &file, uint32_t iso, uint32_t phi, uint32_t etaSign, uint32_t eta, uint32_t et );
      void printJet(ofstream &file, uint32_t phi, uint32_t etaSign, uint32_t eta, uint32_t et );
      void printSum(ofstream &file, uint32_t phi, uint32_t et );
      void printTrack(ofstream &file, float phi, float eta, float et );
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void initializEcalTpgs(edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs,std::vector<TLorentzVector> &allEcalTPGs);
      void initializHcalTpgs(edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs,std::vector<TLorentzVector> &allHcalTPGs,const edm::EventSetup& es);
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      typedef std::vector<TTTrack< Ref_Phase2TrackerDigi_ >> L1TkTrackCollectionType;
      typedef vector<reco::GenParticle> GenParticleCollectionType;

      int run, lumi, event;

      int nev_; // Number of events processed
      bool verbose_;
      std::ofstream logFile_;
      int patternNumber;

      bool compareByPt_tracks (TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j) { 
	return(i.getMomentum().perp() > j.getMomentum().perp()); 
      };
      
      edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalSrc_; 
      edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalSrc_;
      edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
      edm::InputTag L1TrackInputTag;
      edm::InputTag L1TrackPrimaryVertexTag;

      ofstream fin;
      ofstream fout;
      std::string summaryCardOutputFileName_;
      std::string summaryCardInputFileName_;
      double recoPt_;
};


