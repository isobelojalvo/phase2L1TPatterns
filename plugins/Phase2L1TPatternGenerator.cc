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
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//


Phase2L1TPatternGenerator::Phase2L1TPatternGenerator(const edm::ParameterSet& cfg):
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  L1ClustersToken_(consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters")))
{

  L1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
  L1TrackPrimaryVertexTag = cfg.getParameter<edm::InputTag>("L1TrackPrimaryVertexTag");

  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);   

   //now do what ever initialization is needed

  summaryCardOutputFileName_  = cfg.getUntrackedParameter<std::string>("summaryCardOutputFileName");
  summaryCardInputFileName_   = cfg.getUntrackedParameter<std::string>("summaryCardInputFileName");
 
  fout.open(summaryCardOutputFileName_);
  fin.open(summaryCardInputFileName_);

  clustersInputFileName_  = "clustersInput.txt";
  finCluster.open(clustersInputFileName_);

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
  edm::Handle< std::vector<L1CaloCluster> > l1CaloClusters;
  evt.getByToken( L1ClustersToken_, l1CaloClusters);


   std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > l1Tracks;

  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > l1trackHandle;
  evt.getByToken(ttTrackToken_, l1trackHandle);
  
  // L1 Track based primary vertex
  //edm::Handle<L1TkPrimaryVertexCollection> l1PrimaryVertexHandle;
  //evt.getByLabel(L1TrackPrimaryVertexTag, l1PrimaryVertexHandle);
  
  l1Tracks.clear();
  //Find and sort the tracks
  for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
    {
       //edm::Ptr<TTTrack<Ref_PixelDigi_>> ptr(l1trackHandle, track_index);
       edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ptr(l1trackHandle, track_index);
       double pt  = ptr->getMomentum().perp();
       double eta = ptr->getMomentum().eta();

       //only using tracks with eta less than 1.5 and pt greater than 2.5 GeV
       if(abs(eta)<1.5 && pt > 2.5)
	 l1Tracks.push_back(l1trackHandle->at(track_index));       
     }

   std::sort(l1Tracks.begin(), l1Tracks.end(), [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j){return(i.getMomentum().perp() > j.getMomentum().perp());});   

   //decal ecal and hcal tpgs
   edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
   edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  

   std::vector<TLorentzVector> allEcalTPGs;

   if(!evt.getByToken(ecalSrc_, ecalTPGs))
     std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
   else{
     initializEcalTpgs(ecalTPGs , allEcalTPGs);
   }
   
   std::vector<TLorentzVector> allHcalTPGs;

   if(!evt.getByToken(hcalSrc_, hcalTPGs))
     std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
   else{
     initializHcalTpgs(hcalTPGs , allHcalTPGs, iSetup);
   }

   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   
   fout<<"#Event "<<patternNumber<<endl;
   fin <<"#Event "<<patternNumber<<endl;
   patternNumber++;
   std::cout<<"run: "<<run<<" lumi: "<< lumi <<" event: "<< event<<std::endl;

   
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
   
   iWord = 0;
   for(unsigned i = 0; i < l1CaloClusters->size() ; i++){
     L1CaloCluster cluster = l1CaloClusters->at(i);
     //if(cluster.p4().Pt()>0)
     std::cout<<"raw cluster "<<std::hex<<cluster.raw()<<std::endl;
     double pt =  cluster.p4().Pt();
     double eta = cluster.p4().Eta();
     double phi = cluster.p4().Phi();
     std::cout<<"pt "<<pt<<" eta "<<eta<<" phi "<<phi<<std::endl;
     printCluster(finCluster, phi, eta, pt, cluster.raw());
     iWord++;
     if(iWord%10 == 0) finCluster << endl;
     
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


void Phase2L1TPatternGenerator::initializEcalTpgs(edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs,std::vector<TLorentzVector> &allEcalTPGs)
{
  triggerGeometryTools trigTools;
  for (size_t i = 0; i < ecalTPGs->size(); ++i) {
    
    int cal_ieta = (*ecalTPGs)[i].id().ieta();
    int cal_iphi = (*ecalTPGs)[i].id().iphi();
    if(cal_iphi==0)
      std::cout<<"cal_phi is 0"<<std::endl;
    if(cal_ieta<-28)
      continue;
    if(cal_ieta>28)
      continue;
    int ieta = trigTools.TPGEtaRange(cal_ieta);
    short zside = (*ecalTPGs)[i].id().zside();
    // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
    // TPG ieta ideal goes from 0-55.
    double LSB = 0.5;
    double et= (*ecalTPGs)[i].compressedEt()*LSB;
    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
    }

    float eta = trigTools.getRecoEta(ieta, zside);
    float phi = trigTools.getRecoPhi(cal_iphi);    
    //if(et>0)
    //std::cout<<"et "<<et<<std::endl;
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    //if(et>5)
    //std::cout<<"Event Display tpg ecal pt() "<<temp.Pt()<< " eta " <<eta << " phi "<< phi <<std::endl;
    allEcalTPGs.push_back(temp);
  }
}  

void Phase2L1TPatternGenerator::initializHcalTpgs(edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs,std::vector<TLorentzVector> &allHcalTPGs,const edm::EventSetup& es){

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  triggerGeometryTools trigTools;  
  for (size_t i = 0; i < hcalTPGs->size(); ++i) {
    HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    if(cal_ieta>28)continue; 
    if(cal_ieta<-28)continue; 
    int ieta = trigTools.TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside = tpg.id().zside();
    double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //if(et>0)
    //std::cout<<"HCAL ET "<<et<<std::endl;
    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
    }

    float eta = trigTools.getRecoEta(ieta, zside);
    float phi = trigTools.getRecoPhi(cal_iphi);    
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHcalTPGs.push_back(temp);
  }
  
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Phase2L1TPatternGenerator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Phase2L1TPatternGenerator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Phase2L1TPatternGenerator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Phase2L1TPatternGenerator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

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
