// -*- C++ -*-
//
// Package:    Phase2L1TPatternGenerator
// Class:      Phase2L1TPatternGenerator
// 
/**\class Phase2L1TPatternGenerator Phase2L1TPatternGenerator.cc L1Trigger/Phase2L1TPatternGenerator/plugins/Phase2L1TPatternGenerator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Thu, 06 Apr 2017 12:53:23 GMT
// $Id$
//
//


#include "L1Trigger/phase2L1TPatterns/interface/Phase2L1TPatternGenerator.h"
#include <fstream>

#define max_n_tracks 50

using std::endl;
using std::cout;
using std::vector;
using std::setfill;
using std::setw;

Phase2L1TPatternGenerator::Phase2L1TPatternGenerator(const edm::ParameterSet& cfg):
  ecalTPGBToken_(consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalTPGsBarrel"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalTPGs"))),
  L1ClustersToken_(consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters")))
{

  L1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
  L1TrackPrimaryVertexTag = cfg.getParameter<edm::InputTag>("L1TrackPrimaryVertexTag");
  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);   

  summaryCardOutputFileName_  = cfg.getUntrackedParameter<std::string>("summaryCardOutputFileName");
  summaryCardInputFileName_   = cfg.getUntrackedParameter<std::string>("summaryCardInputFileName");
 
  fout.open(summaryCardOutputFileName_);
  fin.open(summaryCardInputFileName_);

  clustersInputFileName_  = "clustersInput.txt";
  finCluster.open(clustersInputFileName_);

  ecalCrystalInputFileName_  = "inputEcalCrystals.txt";
  finCrystal.open(ecalCrystalInputFileName_);

  hcalTPGInputFileName_  = "inputHcalTPGs.txt";
  finHCALTPG.open(hcalTPGInputFileName_);

  std::cout<<"beginning job"<<endl;
  /*
  fout<<"#EG nonIso: 8 highest"<<endl;
  fout<<"#EG Iso: 8 highest"<<endl;
  fout<<"#Tau nonIso: 8 highest"<<endl;
  fout<<"#Tau Iso: 8 highest"<<endl;
  fout<<"#Jets Central: 8 highest"<<endl;
  fout<<"#Jets Highest: 8 highest"<<endl;
  fout<<"#MET "<<endl;
  fout<<"#ET "<<endl;
  fout<<"#MHT "<<endl;
  fout<<"#HT "<<endl;
  fout<<"#PUMBin, pileup"<<endl;
  */
  fout<<"#50 Tracks"<<endl;
  fout<<"#output all TPGs"<<endl;
  fout<<endl;

  patternNumber = 0;
}


Phase2L1TPatternGenerator::~Phase2L1TPatternGenerator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Phase2L1TPatternGenerator::analyze(const edm::Event& evt, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //////////////////// Set up the ECAL Crystals /////////////////////
  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  evt.getByToken( ecalTPGBToken_, ecaltpgCollection);
  vector<ecalCrystal_t> ecalCrystals;

  // Get the ECAL Crystals as ecalCrystal_t
  getEcalCrystals(ecaltpgCollection, ecalCrystals);

  ///////////////////// Set Up the HCAL Crystals /////////////////////
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcalTPGCollection;
  vector<TLorentzVector> hcalTPGs;
  
  if(!evt.getByToken(hcalSrc_, hcalTPGCollection))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;

  edm::ESHandle<L1CaloHcalScale> hcalScale;
  iSetup.get<L1CaloHcalScaleRcd>().get(hcalScale);
  getHcalTPGs(hcalTPGCollection, hcalScale, hcalTPGs);

  ///////////////////// Set Up the Calo Clusters /////////////////////
  edm::Handle< std::vector<L1CaloCluster> > l1CaloClusters;
  evt.getByToken( L1ClustersToken_, l1CaloClusters);

   std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > l1Tracks;

   ///////////////////// Set Up the Tracks /////////////////////
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > l1trackHandle;
  evt.getByToken(ttTrackToken_, l1trackHandle);
  l1Tracks.clear();
  //Find and sort the tracks
  for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
    {
       edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ptr(l1trackHandle, track_index);
       double pt  = ptr->getMomentum().perp();
       double eta = ptr->getMomentum().eta();


       //only using tracks with eta less than 1.5 and pt greater than 2.5 GeV
       if(abs(eta)<1.5 && pt > 2.5)
	 l1Tracks.push_back(l1trackHandle->at(track_index));       
     }

   std::sort(l1Tracks.begin(), 
	     l1Tracks.end(), 
	     [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j)
	     {return(i.getMomentum().perp() > j.getMomentum().perp());});   


   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   
   fout<<"#Event "<<patternNumber<<endl;
   fin <<"#Event "<<patternNumber<<endl;

   patternNumber++;
   std::cout<<"run: "<<run<<" lumi: "<< lumi <<" event: "<< event<<std::endl;

   ///////////////////// Print Tracks /////////////////////
   int iWord = 0;
   //create high pt track distribution
   if(l1Tracks.size()>0){
     for(unsigned int i = 0; i < max_n_tracks && i < l1Tracks.size(); i++){
       TTTrack< Ref_Phase2TrackerDigi_ > l1Track = l1Tracks.at(i);
       double pt = l1Track.getMomentum().perp();
       double eta = l1Track.getMomentum().eta();
       double phi = l1Track.getMomentum().phi();
       printTrack(fin, phi, eta, pt);
       iWord++;
       if(iWord%10 == 0) fin << endl;
     }
   }

   if(l1Tracks.size()<max_n_tracks){
     for(int i = l1Tracks.size(); i < max_n_tracks ; i++){
       printTrack(fin,0,0,0);
     }
   }
   
   ///////////////////// Print Clusters /////////////////////
   iWord = 0;
   for(unsigned i = 0; i < l1CaloClusters->size() ; i++){
     L1CaloCluster cluster = l1CaloClusters->at(i);
     //if(cluster.p4().Pt()>0)
     //std::cout<<"raw cluster "<<std::hex<<cluster.raw()<<std::endl;

     double pt =  cluster.p4().Pt();
     double eta = cluster.p4().Eta();
     double phi = cluster.p4().Phi();
     //std::cout<<"pt "<<pt<<" eta "<<eta<<" phi "<<phi<<std::endl;
     printCluster(finCluster, phi, eta, pt, cluster.raw());
     iWord++;
     if(iWord%10 == 0) finCluster << endl;
   }
   
   ///////////////////// Print Crystals /////////////////////
   iWord = 0;
   for( int iCrystalEta = -100; iCrystalEta < 101; iCrystalEta++){
     for(unsigned int iCrystalPhi = 0; iCrystalPhi < 361; iCrystalPhi++){

       if(iCrystalEta == 0)
	 continue;

       ecalCrystal_t foundCrystal;

       double pt = findEcalCrystal(iCrystalEta, iCrystalPhi, ecalCrystals, foundCrystal);
       //double pt =  foundCrystal.p4.Pt();
       double eta = foundCrystal.p4.Eta();
       double phi = foundCrystal.p4.Phi();

       printECALCrystal(finCrystal, phi, eta, pt);

       iWord++;
       if(iWord%10 == 0) finCrystal << endl;

     }
   }

   ///////////////////// Print HCAL TPGs /////////////////////

  for(int tEta = -20; tEta < 20; tEta++){
    //skip tower eta == 0
    if(tEta == 0)
      continue;

    //loop over all towers in phi
    for(int tPhi = 0; tPhi < 72; tPhi++){
      float pt = 0;
      //tRecoEta is from -1.74 to 1.74 add 0.0435 to get the center of the bin
      float eta = tEta * 0.087 + 0.087/2;      
      //tRecoEta is from 0 to 2pi (bin is 2pi/72) add 0.0436 to get the center of the bin
      float phi = tPhi * 0.0871380 + 0.0871380/2;
      
      //Find matching HcalTPG
      for(auto hcalTPG : hcalTPGs){
	TLorentzVector p4; 
	p4.SetPtEtaPhiE(5,eta,phi,5);
	if(hcalTPG.DeltaR( p4 )< 0.08727/2 ){
	  pt = hcalTPG.Pt();
	  break;
	}
      }
      printHCALTPG(finHCALTPG, phi, eta, pt);
       iWord++;
       if(iWord%10 == 0) finHCALTPG << endl;

    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
Phase2L1TPatternGenerator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Phase2L1TPatternGenerator::endJob() 
{

  fout.close();
  fin.close();
  finCluster.close();
  finHCALTPG.close();
  finCrystal.close();
}

/*
 * The iEta as defined by ECAL are a 10 crystal offset 
 * from the hcal tpgs therefore, search for crystal + 10
 */
float Phase2L1TPatternGenerator::findEcalCrystal(int cEta, int cPhi, 
						vector<ecalCrystal_t> ecalCrystals, 
						ecalCrystal_t &foundCrystal){
  float crystalET = 0;
  for(auto ecalCrystal : ecalCrystals){

    int findEcalIEta = cEta;
    int findEcalIPhi = cPhi + 10; // 10 crystal offset from what is provided as iPhi input to match to hcal

    if(ecalCrystal.iEta == findEcalIEta && 
       ecalCrystal.iPhi == findEcalIPhi){
      crystalET    = ecalCrystal.p4.Pt();
      foundCrystal = ecalCrystal;
      std::cout<<"found ecal crystal "<<crystalET<<std::endl;
    }

  }

  return crystalET;
}


/*
 * 28 bit number
 * |iphi
 * |         |ieta sign
 * |         | 
 * |         ||ieta
 * |         ||         |et
 * 0000 0000 0000 0000 0000 0000 0000
 */

void Phase2L1TPatternGenerator::printTrack(ofstream &file, float phi, float eta, float et ){
  triggerGeometryTools trigTools;
  uint32_t et_int  = round(et*10); 
  uint32_t eta_uint = trigTools.getCrystalIEta(eta); 
  uint32_t phi_uint = trigTools.getCrystalIPhi(phi); 
  
  if(et>200)
    et_int = 0x7FF;
  if(et > 0x7FF){
    std::cout<<"you are trying to print a track with pt value greater than 0xFFF, something is wrong."<<std::endl;
  }

  uint32_t etaSign = 0;
  if(eta/abs(eta)>0)
    etaSign = 1;
  else
    etaSign = 0;

  uint32_t printNumber = 0;
  printNumber += et_int&0x7FF;
  printNumber += ((eta_uint&0x7F)<<11);
  printNumber += (etaSign<<19);
  printNumber += ((phi_uint&0xFF)<<20);
  file<<std::hex<<setfill('0') << setw(8)<<printNumber<<" ";
  
}

/*
 * ECAL Crystal
 *   |et
 *  0000 0000 0000
 */

void Phase2L1TPatternGenerator::printECALCrystal(ofstream &file, float phi, float eta, float et ){
  triggerGeometryTools trigTools;
  uint32_t et_int  = round(et*10); 
  //uint32_t eta_uint = trigTools.getCrystalIEta(eta); 
  //uint32_t phi_uint = trigTools.getCrystalIPhi(phi); 
  
  if(et>200)
    et_int = 0x7FF;
  if(et > 0x7FF){
    std::cout<<"you are trying to print a track with pt value greater than 0xFFF, something is wrong."<<std::endl;
  }

  uint32_t printNumber = 0;
  printNumber += et_int & 0x7FF;
  file<<std::hex<<setfill('0') << setw(4)<<printNumber<<" ";
  
}


/*
 * ECAL Crystal
 *   |et
 *  0000 0000 0000
 */

void Phase2L1TPatternGenerator::printHCALTPG(ofstream &file, float phi, float eta, float et ){
  triggerGeometryTools trigTools;
  uint32_t et_int  = round(et*10); 
  //uint32_t eta_uint = trigTools.getCrystalIEta(eta); 
  //uint32_t phi_uint = trigTools.getCrystalIPhi(phi); 
  
  if(et>200)
    et_int = 0x7FF;
  if(et > 0x7FF){
    std::cout<<"you are trying to print a track with pt value greater than 0x7FF, something is wrong."<<std::endl;
  }

  uint32_t printNumber = 0;
  printNumber += et_int & 0x7FF;
  file<<std::hex<<setfill('0') << setw(4)<<printNumber<<" ";
  
}


/*
 * 28 bit number
 * |iphi
 * |         |ieta sign
 * |         | 
 * |         ||ieta
 * |         ||         |et
 * 0000 0000 0000 0000 0000 0000 0000
 */

void Phase2L1TPatternGenerator::printCluster(ofstream &file, float phi, float eta, float et, uint32_t raw ){
  triggerGeometryTools trigTools;
  uint32_t et_int  = round(et*10); 
  uint32_t eta_uint = trigTools.getCrystalIEta(eta); 
  uint32_t phi_uint = trigTools.getCrystalIPhi(phi); 
  
  if(et>200)
    et_int = 0x7FF;
  if(et > 0x7FF){
    std::cout<<"you are trying to print a track with pt value greater than 0xFFF, something is wrong."<<std::endl;
  }

  uint32_t etaSign = 0;
  if(eta/abs(eta)>0)
    etaSign = 1;
  else
    etaSign = 0;

  uint32_t printNumber = 0;
  printNumber += et_int&0x7FF;
  printNumber += ((eta_uint&0x7F)<<11);
  printNumber += (etaSign<<19);
  printNumber += ((phi_uint&0xFF)<<20);
  file<<std::hex<<setfill('0') << setw(8)<<printNumber<<" ";
  
}

void Phase2L1TPatternGenerator::printEGTau(ofstream &file, uint32_t iso, uint32_t phi, uint32_t etaSign, uint32_t eta, uint32_t et ){
  uint32_t printNumber = 0;
  printNumber += et&0x1FF;
  printNumber += ((eta&0xFF)<<9);
  printNumber += (etaSign<<16);
  printNumber += ((phi&0xFF)<<17);
  printNumber += ((iso&0x1)<<25);
  file<<std::hex<<setfill('0') << setw(8)<<printNumber<<" ";
}

void Phase2L1TPatternGenerator::printJet(ofstream &file, uint32_t phi, uint32_t etaSign, uint32_t eta, uint32_t et ){
  uint32_t printNumber = 0;
  printNumber += et&0x3FF;
  printNumber += ((eta&0xFF)<<11);
  printNumber += (etaSign<<18);
  printNumber += ((phi&0xFF)<<19);
  file<<std::hex<<setfill('0') << setw(8)<<printNumber<<" ";
}

void Phase2L1TPatternGenerator::printSum(ofstream &file, uint32_t phi, uint32_t et ){
  uint32_t printNumber = 0;
  printNumber += et&0x7FF;
  printNumber += ((phi&0xFF)<<9);
  file<<std::hex<<setfill('0') << setw(8)<<printNumber<<" ";
}


void Phase2L1TPatternGenerator::getHcalTPGs( edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcaltpgCollection, 
				   edm::ESHandle<L1CaloHcalScale> &hcalScale, 
				   vector<TLorentzVector> &allHcalTPGs){

  for (size_t i = 0; i < hcaltpgCollection->size(); ++i) {

    HcalTriggerPrimitiveDigi tpg = (*hcaltpgCollection)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    if(cal_ieta> 28)continue; 
    if(cal_ieta<-28)continue; 

    int ieta      = TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside   = tpg.id().zside();
    double et     = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //DEBUG STATEMENT
    //if(et>0)
    //std::cout<<"HCAL ET "<<et<<std::endl;

    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
    }

    triggerGeometryTools tool;
    float eta = tool.getRecoEta(ieta, zside);
    float phi = tool.getRecoPhi(cal_iphi);    

    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);

    allHcalTPGs.push_back(temp);
  }
  
}

void Phase2L1TPatternGenerator::getEcalCrystals(edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs, vector<ecalCrystal_t> &ecalCrystals)
{
  
  for(auto& tpg : *ecalTPGs.product())
    {
      if(tpg.encodedEt() > 0) 
	{

	  GlobalVector position;
	  auto cell = ebGeometry->getGeometry(tpg.id());

	  float et = tpg.encodedEt()/8.;

	  if(et<0.5) continue;
	  //float energy = et / sin(position.theta());
	  float eta = cell->getPosition().eta();
	  float phi = cell->getPosition().phi();
	  EBDetId detID = tpg.id();
	  float iEta = detID.ieta();
	  float iPhi = detID.iphi();

	  //DEBUG STATMENT
	  /*
	  if(et>0.5){
	    std::cout<<"ET "<<et<<std::endl;
	    std::cout<<"eta  "<< eta<< " phi  "<< phi<<std::endl;
	    std::cout<<"iEta"<<iEta<<" iphi "<<iPhi<<std::endl;
	    }*/
	  ecalCrystal_t tempCrystal;
	  tempCrystal.p4.SetPtEtaPhiE(et, eta, phi,et);
	  tempCrystal.iEta = iEta;
	  tempCrystal.iPhi = iPhi;
	  tempCrystal.id = detID;

	  ecalCrystals.push_back(tempCrystal);
	}
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Phase2L1TPatternGenerator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1TPatternGenerator);
