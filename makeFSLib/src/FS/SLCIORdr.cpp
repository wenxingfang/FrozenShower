#include "SLCIORdr.h"

#include "lcio.h"  //LCIO
#include "LCIOSTLTypes.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCIO.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "IO/LCReader.h"
#include "IMPL/MCParticleImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "UTIL/CellIDDecoder.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#include <map>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include "my_utils.h"

#ifndef PI
#define PI acos(-1)
#endif
#ifndef DEBUG
#define DEBUG true
#endif
using namespace lcio;
using namespace IMPL;
using namespace std;

// add HCALBarrel hit //

SLCIORdr::SLCIORdr(string name, string output, bool is_gun, int pid){


m_slcio_rdr = IOIMPL::LCFactory::getInstance()->createLCReader();
m_slcio_rdr->open(name.c_str());
m_processed_event=0;
m_is_gun = is_gun;
m_pid    = pid   ;
file_out = new TFile(output.c_str(),"RECREATE"); 
tree_out = new TTree("evt","test tree");
m_HitCm_x   = 0;
m_HitCm_y   = 0;
m_HitCm_z   = 0;
m_HitEn_tot = 0;

//tree_out->Branch("evt_id",&m_processed_event,"evt_id/I");
//tree_out->Branch("m_HitCm_x"  , &m_HitCm_x   ,"m_HitCm_x/D"     );
//tree_out->Branch("m_HitCm_y"  , &m_HitCm_y   ,"m_HitCm_y/D"     );
//tree_out->Branch("m_HitCm_z"  , &m_HitCm_z   ,"m_HitCm_z/D"     );
//tree_out->Branch("m_HitEn_tot", &m_HitEn_tot ,"m_HitEn_tot/D"   );
//tree_out->Branch("m_HitFirst_x"     , &m_HitFirst_x      );
//tree_out->Branch("m_HitFirst_y"     , &m_HitFirst_y      );
//tree_out->Branch("m_HitFirst_z"     , &m_HitFirst_z      );
//tree_out->Branch("m_phi_rotated"    , &m_phi_rotated     );
//tree_out->Branch("m_HitFirst_vtheta", &m_HitFirst_vtheta );
//tree_out->Branch("m_HitFirst_vphi"  , &m_HitFirst_vphi   );
tree_out->Branch("m_mc_vertexX", &m_mc_vertexX  );
tree_out->Branch("m_mc_vertexY", &m_mc_vertexY  );
tree_out->Branch("m_mc_vertexZ", &m_mc_vertexZ  );
tree_out->Branch("m_mc_Px", &m_mc_Px       );
tree_out->Branch("m_mc_Py", &m_mc_Py       );
tree_out->Branch("m_mc_Pz", &m_mc_Pz       );
tree_out->Branch("m_mc_M", &m_mc_M        );
tree_out->Branch("m_Hit_x", &m_Hit_x       );
tree_out->Branch("m_Hit_y", &m_Hit_y       );
tree_out->Branch("m_Hit_z", &m_Hit_z       );
tree_out->Branch("m_Hit_E", &m_Hit_E       ); 
tree_out->Branch("m_HcalHit_x", &m_HcalHit_x       );
tree_out->Branch("m_HcalHit_y", &m_HcalHit_y       );
tree_out->Branch("m_HcalHit_z", &m_HcalHit_z       );
tree_out->Branch("m_HcalHit_E", &m_HcalHit_E       ); 
tree_out->Branch("m_mc_Pdg", &m_mc_Pdg      );
//tree_out->Branch("m_mc_Charge", &m_mc_Charge      );
//tree_out->Branch("m_mc_genStatus", &m_mc_genStatus); 
//tree_out->Branch("m_mc_simStatus", &m_mc_simStatus); 
//tree_out->Branch("m_mc_np", &m_mc_np       ); 
//tree_out->Branch("m_mc_nd", &m_mc_nd       ); 
//tree_out->Branch("m_Hit_ID0", &m_Hit_ID0      );
//tree_out->Branch("m_Hit_ID1", &m_Hit_ID1      );
tree_out->Branch("m_ID_S", &m_ID_S      );
tree_out->Branch("m_ID_M", &m_ID_M      );
tree_out->Branch("m_ID_I", &m_ID_I      );
tree_out->Branch("m_ID_J", &m_ID_J      );
tree_out->Branch("m_ID_K", &m_ID_K      );
//tree_out->Branch("m_mc_pHitx"       , &m_mc_pHitx       );
//tree_out->Branch("m_mc_pHity"       , &m_mc_pHity       );
//tree_out->Branch("m_mc_pHitz"       , &m_mc_pHitz       );
//tree_out->Branch("m_mc_pHit_theta"  , &m_mc_pHit_theta  );
//tree_out->Branch("m_mc_pHit_phi"    , &m_mc_pHit_phi    );
//tree_out->Branch("m_mc_pHit_rotated", &m_mc_pHit_rotated);
//tree_out->Branch("m_Hits", &m_Hits      ); 
//tree_out->Branch("m_HcalHits", &m_HcalHits      ); 
//tree_out->Branch("m_mc_pHit_dz", &m_mc_pHit_dz      ); 
//tree_out->Branch("m_mc_pHit_dy", &m_mc_pHit_dy      ); 


std::cout<<"pi="<< PI << std::endl;

}

SLCIORdr::~SLCIORdr(){
delete m_slcio_rdr;
}

bool SLCIORdr::mutate(){

    clear();
    EVENT::LCEvent *lcEvent = m_slcio_rdr->readNextEvent(LCIO::UPDATE);
    LCCollection *Col = NULL;
    //LCCollectionVec *ColVec = NULL;
    LCCollection *ColVec = NULL;

    if(lcEvent) {  //cout<<"det name="<< lcEvent->getDetectorName() <<endl;
                   const std::vector<std::string>* colVec = lcEvent->getCollectionNames();
                   /*
                   int skip = 0;
                   for(unsigned int i=0; i< colVec->size(); i++)
                   {
                       //if(colVec->at(i)=="ECALBarrel" || colVec->at(i)=="HCALBarrel"){skip = skip+1 ;}
                       if(colVec->at(i)=="EcalBarrelSiliconCollection"){skip = skip+1 ;}
                       //std::cout<<m_processed_event<<";"<<colVec->at(i)<<std::endl;
                   }
                   if (skip != 1) return true;
                   */
                   for(unsigned int i=0; i< colVec->size(); i++){
                       //cout<<"col name="<< colVec->at(i) <<endl;
                       if(colVec->at(i)=="EcalBarrelSiliconCollection"){// consider only barrel now
                       //if(colVec->at(i)=="EcalEndcapSiliconCollection"){// consider only endcap now
                       //if(colVec->at(i)=="EcalBarrelSiliconPreShowerCollection"){// consider only barrel now
                       ColVec  = lcEvent->getCollection(colVec->at(i));
                       CellIDDecoder<SimCalorimeterHit> IdDecoder(ColVec);
                       getSimHitInfo(IdDecoder, ColVec, m_Hit_ID0, m_Hit_ID1, m_Hit_x, m_Hit_y, m_Hit_z, m_Hit_E , m_HitCm_x, m_HitCm_y, m_HitCm_z, m_HitEn_tot, m_ID_S, m_ID_M, m_ID_I, m_ID_J, m_ID_K );
                       }
                       else if(colVec->at(i)=="HcalBarrelCollection"){// consider only barrel now
                       //else if(colVec->at(i)=="HcalEndcapCollection"){// consider only endcap now
                       ColVec  = lcEvent->getCollection(colVec->at(i));
                       getSimHcalHitInfo(ColVec, m_HcalHit_x, m_HcalHit_y, m_HcalHit_z, m_HcalHit_E);
                       }
                       else if (colVec->at(i)=="MCParticle")
                       {
                           Col = lcEvent->getCollection("MCParticle");
                           getMCinfo0(Col, m_mc_Pdg, m_mc_Charge,m_mc_genStatus, m_mc_simStatus, m_mc_vertexX, m_mc_vertexY, m_mc_vertexZ, m_mc_Px, m_mc_Py, m_mc_Pz, m_mc_M);
                       }// if mc particle
                   }
                  }
    else {cout<<"end file, total event="<<m_processed_event<<"\n"; return false;}
    m_processed_event ++;
    if(m_mc_M.size() == 1)tree_out->Fill();
    if(m_processed_event%10000==0)cout << "Done for event "<<m_processed_event<< endl;
    return true;
}

bool SLCIORdr::getMCinfo(LCCollection* lcCol, std::vector<int>& mc_Pdg, std::vector<double>& mc_Charge, std::vector<int>& mc_genStatus, std::vector<int>& mc_simStatus, std::vector<double>& mc_vertexX, std::vector<double>& mc_vertexY, std::vector<double>& mc_vertexZ, std::vector<double>& mc_Px, std::vector<double>& mc_Py, std::vector<double>& mc_Pz, std::vector<double>& mc_M , std::vector<int>& mc_np, std::vector<int>& mc_nd, std::vector<double>& mc_pHitx, std::vector<double>& mc_pHity, std::vector<double>& mc_pHitz, std::vector<double>& mc_pHit_theta, std::vector<double>& mc_pHit_phi, std::vector<double>& mc_pHit_rotated, LCCollection* SimHitCol, CellIDDecoder<SimCalorimeterHit> & SimHit_idDecoder, LCCollection* DigiHitCol, CellIDDecoder<CalorimeterHit> & DigiHit_idDecoder, std::vector< std::vector<double> >& Hits, std::vector<double>& mc_pHit_dz, std::vector<double>& mc_pHit_dy, LCCollection* DigiHcalCol, std::vector< std::vector<double> >& HcalHits){
      
    if(lcCol){
        //if(DEBUG)std::cout<<"Hi:"<<std::endl;
	int NHEP = lcCol->getNumberOfElements();
	for( int IHEP=0; IHEP<NHEP; IHEP++ ){
	  EVENT::MCParticle* in = dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP));
          //if(fabs(in->getPDG())!=m_pid ) continue;
          if(in->getPDG() != m_pid ) continue;
	  EVENT::MCParticleVec parents = (dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP)))->getParents();
	  EVENT::MCParticleVec daughters = (dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP)))->getDaughters();
	  int np = parents.size();
	  int nd = daughters.size();
          if(m_is_gun && np!=0) continue;
          //std::cout<<"check1:"<<in->getPDG()<<<<","<<m_is_gun<<","<<np<<std::endl; 
          //float Distance = 1850; //mm
          float Distance = 1840; //mm, remove some back scatter particles
          float BField = 3.0; //T
	  //if (np != 0 ) continue; // only select original particle 
	  //if (in->isCreatedInSimulation() == true) continue; // only select gen particle 
	  //if (in->getGeneratorStatus() != 1) continue; // only select stable gen particle 
	  //if (fabs(in->getPDG()) != 11) continue; // only select e+- 
          bool vertex_beforeEcal   = beforeECAL(in->getVertex()[0]  , in->getVertex()[1]  , Distance);
          bool endpoint_beforeEcal = beforeECAL(in->getEndpoint()[0], in->getEndpoint()[1], Distance);
          if(vertex_beforeEcal!=true || endpoint_beforeEcal!=false) continue;
          //if(in->getCharge()==0) continue;           
          //if(fabs(in->getPDG())!=11 ) continue;
          float pHitx=1850;// mm the x place for ECAL
          float pHity=0;
          float pHitz=0;
          float pHit_theta=0;
          float pHit_phi=0;
          float pHit_rotated=0;
          int getHit = getHitPoint(in->getCharge(), in->getVertex()[0], in->getVertex()[1], in->getVertex()[2], in->getMomentum()[0], in->getMomentum()[1], in->getMomentum()[2], BField, in->getEndpoint()[0], in->getEndpoint()[1], pHitx, pHity, pHitz, pHit_theta, pHit_phi, pHit_rotated);
          if(fabs(pHity) > 600) continue;//remove high Y region now
          if(fabs(pHitz) > 2300)continue;//remove endcap now
          //std::cout<<"check2:"<<"pHitx="<<pHitx<<",pHity="<<pHity<<",pHitz="<<pHitz<<std::endl; 
	  //cout << "DEBUG: PDG= " << in->getPDG() <<",pHitX="<<pHitx<<",pHitY="<<pHity<<",pHitZ="<<pHitz<<",rotated="<<pHit_rotated<<",EndX="<<in->getEndpoint()[0]<<",EndY="<<in->getEndpoint()[1]<<",EndZ="<<in->getEndpoint()[2]<< ", Mx= " << in->getMomentum()[0]<<", My= "<<in->getMomentum()[1]<<",Mz= " <<in->getMomentum()[2]<< "\n ";
          float real_z = pHitz ;
          float real_r = sqrt(pHitx*pHitx + pHity*pHity);
          float real_phi = getPhi(pHitx, pHity)+ pHit_rotated; 
          if(real_phi>=360) real_phi = real_phi - 360 ;
          float real_x = real_r*cos(real_phi*PI/180); 
          float real_y = real_r*sin(real_phi*PI/180); 
          std::vector<double>  tmp_Hits;
          std::vector<double>  tmp_HitsHcal;
          float cell_x = 0 ;
          float cell_y = 0 ;
          float cell_z = 0 ;
          bool access = getHits(SimHitCol, SimHit_idDecoder, DigiHitCol, DigiHit_idDecoder, real_x, real_y, real_z, tmp_Hits, cell_x, cell_y, cell_z);
          const float r_min = 2050;
          const float r_max = 3150;
          const int r_bin = 55;//20 mm each bin
          const float z_half = 600;//+-
          const int z_bin = 20;//30 mm each bin
          const float phi_half = 10;//+-
          const int phi_bin = 20;//0.5  each bin
          
          bool access1 = getHcalHits(DigiHcalCol, real_x, real_y, real_z, tmp_HitsHcal, r_min, r_max, r_bin, z_half, z_bin, phi_half, phi_bin);
          double tmp_dz = pHitz - cell_z;
          float cell_phi = getPhi(cell_x, cell_y) - pHit_rotated;
          double tmp_dy = pHity - sqrt(cell_x*cell_x + cell_y*cell_y)*sin(cell_phi*PI/180);
          if(sqrt(tmp_dz*tmp_dz + tmp_dy*tmp_dy)>5*sqrt(2)) {std::cout<<"too large"<<std::endl;continue;}
          //std::cout<<"tmp_Hits size="<<tmp_Hits.size()<<std::endl;
          mc_Pdg      .push_back(in->getPDG());
          mc_Charge   .push_back(in->getCharge());
	  mc_genStatus.push_back(in->getGeneratorStatus());
	  mc_simStatus.push_back(in->getSimulatorStatus());
	  mc_vertexX  .push_back(in->getVertex()[0]);
	  mc_vertexY  .push_back(in->getVertex()[1]);
	  mc_vertexZ  .push_back(in->getVertex()[2]);
	  mc_Px       .push_back(in->getMomentum()[0]);
	  mc_Py       .push_back(in->getMomentum()[1]);
	  mc_Pz       .push_back(in->getMomentum()[2]);
	  mc_M        .push_back(in->getMass());
	  mc_np       .push_back(np);
	  mc_nd       .push_back(nd);


	  mc_pHitx       .push_back(pHitx);
	  mc_pHity       .push_back(pHity);
	  mc_pHitz       .push_back(pHitz);
          if (pHit_theta < 0) pHit_theta = 180 + pHit_theta ;
	  mc_pHit_theta  .push_back(pHit_theta);
	  mc_pHit_phi    .push_back(pHit_phi);
	  mc_pHit_rotated.push_back(pHit_rotated);

	  Hits.push_back(tmp_Hits);
	  HcalHits.push_back(tmp_HitsHcal);
          mc_pHit_dz.push_back(tmp_dz);
          mc_pHit_dy.push_back(tmp_dy);
          //cout<<"pHit_rotated="<<pHit_rotated<<",hit_x="<<pHitx<<",hit_y="<<pHity<<",hit_z="<<pHitz<<",real_x="<<real_x<<",real_y="<<real_y<<",real_z="<<real_z<<"cell_x="<<cell_x<<",cell_y="<<cell_y<<",cell_z="<<cell_z<<endl;
          //cout<<"tmp_dz="<<tmp_dz<<",tmp_dy="<<tmp_dy<<endl;
          
	}
      }
      else {
	  //cout << "Debug: no MCParticle Collection is read!" << endl;
          return false;
      }
      return true;
}

void SLCIORdr::getMCinfo0(LCCollection* lcCol, std::vector<int>& mc_Pdg, std::vector<double>& mc_Charge, std::vector<int>& mc_genStatus, std::vector<int>& mc_simStatus, std::vector<double>& mc_vertexX, std::vector<double>& mc_vertexY, std::vector<double>& mc_vertexZ, std::vector<double>& mc_Px, std::vector<double>& mc_Py, std::vector<double>& mc_Pz, std::vector<double>& mc_M){
      
    if(lcCol){
        //if(DEBUG)std::cout<<"Hi:"<<std::endl;
	int NHEP = lcCol->getNumberOfElements();
	for( int IHEP=0; IHEP<NHEP; IHEP++ ){
	  EVENT::MCParticle* in = dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP));
          //if(fabs(in->getPDG())!=m_pid ) continue;
          //if(in->getPDG() != m_pid ) continue;
	  EVENT::MCParticleVec parents = (dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP)))->getParents();
	  EVENT::MCParticleVec daughters = (dynamic_cast<EVENT::MCParticle*>(lcCol->getElementAt(IHEP)))->getDaughters();
	  int np = parents.size();
	  int nd = daughters.size();
          if(m_is_gun && np!=0) continue;
          /*
          //std::cout<<"check1:"<<in->getPDG()<<<<","<<m_is_gun<<","<<np<<std::endl; 
          //float Distance = 1850; //mm
          float Distance = 1840; //mm, remove some back scatter particles
          float BField = 3.0; //T
	  //if (np != 0 ) continue; // only select original particle 
	  //if (in->isCreatedInSimulation() == true) continue; // only select gen particle 
	  //if (in->getGeneratorStatus() != 1) continue; // only select stable gen particle 
	  //if (fabs(in->getPDG()) != 11) continue; // only select e+- 
          bool vertex_beforeEcal   = beforeECAL(in->getVertex()[0]  , in->getVertex()[1]  , Distance);
          bool endpoint_beforeEcal = beforeECAL(in->getEndpoint()[0], in->getEndpoint()[1], Distance);
          if(vertex_beforeEcal!=true || endpoint_beforeEcal!=false) continue;
          //if(in->getCharge()==0) continue;           
          //if(fabs(in->getPDG())!=11 ) continue;
          float pHitx=1850;// mm the x place for ECAL
          float pHity=0;
          float pHitz=0;
          float pHit_theta=0;
          float pHit_phi=0;
          float pHit_rotated=0;
          int getHit = getHitPoint(in->getCharge(), in->getVertex()[0], in->getVertex()[1], in->getVertex()[2], in->getMomentum()[0], in->getMomentum()[1], in->getMomentum()[2], BField, in->getEndpoint()[0], in->getEndpoint()[1], pHitx, pHity, pHitz, pHit_theta, pHit_phi, pHit_rotated);
          if(fabs(pHity) > 600) continue;//remove high Y region now
          if(fabs(pHitz) > 2300)continue;//remove endcap now
          //std::cout<<"check2:"<<"pHitx="<<pHitx<<",pHity="<<pHity<<",pHitz="<<pHitz<<std::endl; 
	  //cout << "DEBUG: PDG= " << in->getPDG() <<",pHitX="<<pHitx<<",pHitY="<<pHity<<",pHitZ="<<pHitz<<",rotated="<<pHit_rotated<<",EndX="<<in->getEndpoint()[0]<<",EndY="<<in->getEndpoint()[1]<<",EndZ="<<in->getEndpoint()[2]<< ", Mx= " << in->getMomentum()[0]<<", My= "<<in->getMomentum()[1]<<",Mz= " <<in->getMomentum()[2]<< "\n ";
          float real_z = pHitz ;
          float real_r = sqrt(pHitx*pHitx + pHity*pHity);
          float real_phi = getPhi(pHitx, pHity)+ pHit_rotated; 
          if(real_phi>=360) real_phi = real_phi - 360 ;
          float real_x = real_r*cos(real_phi*PI/180); 
          float real_y = real_r*sin(real_phi*PI/180); 
          std::vector<double>  tmp_Hits;
          std::vector<double>  tmp_HitsHcal;
          float cell_x = 0 ;
          float cell_y = 0 ;
          float cell_z = 0 ;
          bool access = getHits(SimHitCol, SimHit_idDecoder, DigiHitCol, DigiHit_idDecoder, real_x, real_y, real_z, tmp_Hits, cell_x, cell_y, cell_z);
          const float r_min = 2050;
          const float r_max = 3150;
          const int r_bin = 55;//20 mm each bin
          const float z_half = 600;//+-
          const int z_bin = 20;//30 mm each bin
          const float phi_half = 10;//+-
          const int phi_bin = 20;//0.5  each bin
          
          bool access1 = getHcalHits(DigiHcalCol, real_x, real_y, real_z, tmp_HitsHcal, r_min, r_max, r_bin, z_half, z_bin, phi_half, phi_bin);
          double tmp_dz = pHitz - cell_z;
          float cell_phi = getPhi(cell_x, cell_y) - pHit_rotated;
          double tmp_dy = pHity - sqrt(cell_x*cell_x + cell_y*cell_y)*sin(cell_phi*PI/180);
          if(sqrt(tmp_dz*tmp_dz + tmp_dy*tmp_dy)>5*sqrt(2)) {std::cout<<"too large"<<std::endl;continue;}
          //std::cout<<"tmp_Hits size="<<tmp_Hits.size()<<std::endl;
          */
          mc_Pdg      .push_back(in->getPDG());
          mc_Charge   .push_back(in->getCharge());
	  mc_genStatus.push_back(in->getGeneratorStatus());
	  mc_simStatus.push_back(in->getSimulatorStatus());
	  mc_vertexX  .push_back(in->getVertex()[0]);
	  mc_vertexY  .push_back(in->getVertex()[1]);
	  mc_vertexZ  .push_back(in->getVertex()[2]);
	  mc_Px       .push_back(in->getMomentum()[0]);
	  mc_Py       .push_back(in->getMomentum()[1]);
	  mc_Pz       .push_back(in->getMomentum()[2]);
	  mc_M        .push_back(in->getMass());
          
	}
      }
}



bool SLCIORdr::getDigiHitInfo(CellIDDecoder<CalorimeterHit> & Hit_idDecoder, LCCollection* Col, std::vector<int>& vec_ID0, std::vector<int>& vec_ID1, std::vector<double>& vec_Hit_x, std::vector<double>& vec_Hit_y, std::vector<double>& vec_Hit_z, std::vector<double>& vec_Hit_E , double& Hit_cm_x, double& Hit_cm_y, double& Hit_cm_z, double& Hit_tot_e, std::vector<int>& vec_ID_S, std::vector<int>& vec_ID_M, std::vector<int>& vec_ID_I, std::vector<int>& vec_ID_J, std::vector<int>& vec_ID_K){
      if(Col)
      { 
         double tot_E = 0;
         double sum_x = 0;
         double sum_y = 0;
         double sum_z = 0;
         int NHit = Col->getNumberOfElements();
         for( int i=0; i<NHit; i++ ){
           IMPL::CalorimeterHitImpl* in = dynamic_cast<IMPL::CalorimeterHitImpl*>(Col->getElementAt(i));
          //cout<<"Hit "<<i <<",CellID0="<<in->getCellID0()<<", CellID1="<<in->getCellID1()<<", E= "<<in->getEnergy()<<" GeV, x="<<in->getPosition()[0]<<",y="<<in->getPosition()[1]<<",z="<<in->getPosition()[2]<<"\n";
          int StaveNum = Hit_idDecoder(in)["S-1"];// from 0 - 7
          int MNum     = Hit_idDecoder(in)["M"  ];// from 1 - 5
          int ZNum     = Hit_idDecoder(in)["J"  ];// from 0 - 89
          int INum     = Hit_idDecoder(in)["I"  ];// from 0 - 169 for layer 0, then decrease for higher layer
          int LayerNum = Hit_idDecoder(in)["K-1"];// from 0 - 28
          
          vec_ID_S.push_back(StaveNum);
          vec_ID_M.push_back(MNum    );
          vec_ID_I.push_back(INum    );
          vec_ID_J.push_back(ZNum    );
          vec_ID_K.push_back(LayerNum);

          vec_ID0.push_back(in->getCellID0());
          vec_ID1.push_back(in->getCellID1());
          vec_Hit_x.push_back(in->getPosition()[0]);
          vec_Hit_y.push_back(in->getPosition()[1]);
          vec_Hit_z.push_back(in->getPosition()[2]);
          vec_Hit_E.push_back(in->getEnergy());
          tot_E = tot_E + in->getEnergy() ;
          sum_x = sum_x + (in->getEnergy())*(in->getPosition()[0]);
          sum_y = sum_y + (in->getEnergy())*(in->getPosition()[1]);
          sum_z = sum_z + (in->getEnergy())*(in->getPosition()[2]);
          }
         Hit_cm_x = Hit_cm_x + sum_x; 
         Hit_cm_y = Hit_cm_y + sum_y; 
         Hit_cm_z = Hit_cm_z + sum_z; 
         Hit_tot_e = Hit_tot_e + tot_E; 
         //std::cout<<"digi tot_e="<<Hit_tot_e<<std::endl;
      }
      else {
	  //cout << "Debug: no "<<Col->getTypeName()<<" Collection is found!" << endl;
          return false;
      }
    return true;

}

bool SLCIORdr::getSimHitInfo(CellIDDecoder<SimCalorimeterHit> & Hit_idDecoder, LCCollection* Col, std::vector<int>& vec_ID0, std::vector<int>& vec_ID1, std::vector<double>& vec_Hit_x, std::vector<double>& vec_Hit_y, std::vector<double>& vec_Hit_z, std::vector<double>& vec_Hit_E , double& Hit_cm_x, double& Hit_cm_y, double& Hit_cm_z, double& Hit_tot_e, std::vector<int>& vec_ID_S, std::vector<int>& vec_ID_M, std::vector<int>& vec_ID_I, std::vector<int>& vec_ID_J, std::vector<int>& vec_ID_K){
      if(Col)
      { 
         double tot_E = 0;
         double sum_x = 0;
         double sum_y = 0;
         double sum_z = 0;
         int NHit = Col->getNumberOfElements();
         for( int i=0; i<NHit; i++ ){
           EVENT::SimCalorimeterHit* in = dynamic_cast<EVENT::SimCalorimeterHit*>(Col->getElementAt(i));
           //IMPL::CalorimeterHitImpl* in = dynamic_cast<IMPL::CalorimeterHitImpl*>(Col->getElementAt(i));
          //cout<<"Hit "<<i <<",CellID0="<<in->getCellID0()<<", CellID1="<<in->getCellID1()<<", E= "<<in->getEnergy()<<" GeV, x="<<in->getPosition()[0]<<",y="<<in->getPosition()[1]<<",z="<<in->getPosition()[2]<<"\n";
          int StaveNum = Hit_idDecoder(in)["S-1"];// from 0 - 7
          int MNum     = Hit_idDecoder(in)["M"  ];// from 1 - 5
          int ZNum     = Hit_idDecoder(in)["J"  ];// from 0 - 89
          int INum     = Hit_idDecoder(in)["I"  ];// from 0 - 169 for layer 0, then decrease for higher layer
          int LayerNum = Hit_idDecoder(in)["K-1"];// from 0 - 28
          
          vec_ID_S.push_back(StaveNum);
          vec_ID_M.push_back(MNum    );
          vec_ID_I.push_back(INum    );
          vec_ID_J.push_back(ZNum    );
          vec_ID_K.push_back(LayerNum);

          vec_ID0.push_back(in->getCellID0());
          vec_ID1.push_back(in->getCellID1());
          vec_Hit_x.push_back(in->getPosition()[0]);
          vec_Hit_y.push_back(in->getPosition()[1]);
          vec_Hit_z.push_back(in->getPosition()[2]);
          vec_Hit_E.push_back(in->getEnergy());
          tot_E = tot_E + in->getEnergy() ;
          sum_x = sum_x + (in->getEnergy())*(in->getPosition()[0]);
          sum_y = sum_y + (in->getEnergy())*(in->getPosition()[1]);
          sum_z = sum_z + (in->getEnergy())*(in->getPosition()[2]);
          }
         Hit_cm_x = Hit_cm_x + sum_x; 
         Hit_cm_y = Hit_cm_y + sum_y; 
         Hit_cm_z = Hit_cm_z + sum_z; 
         Hit_tot_e = Hit_tot_e + tot_E; 
         //std::cout<<"digi tot_e="<<Hit_tot_e<<std::endl;
      }
      else {
	  //cout << "Debug: no "<<Col->getTypeName()<<" Collection is found!" << endl;
          return false;
      }
    return true;

}

bool SLCIORdr::getSimHcalHitInfo(LCCollection* Col, std::vector<double>& vec_Hit_x, std::vector<double>& vec_Hit_y, std::vector<double>& vec_Hit_z, std::vector<double>& vec_Hit_E){
      if(Col)
      { 
         int NHit = Col->getNumberOfElements();
         for( int i=0; i<NHit; i++ ){
           EVENT::SimCalorimeterHit* in = dynamic_cast<EVENT::SimCalorimeterHit*>(Col->getElementAt(i));
          vec_Hit_x.push_back(in->getPosition()[0]);
          vec_Hit_y.push_back(in->getPosition()[1]);
          vec_Hit_z.push_back(in->getPosition()[2]);
          vec_Hit_E.push_back(in->getEnergy());
          }
      }
      else {
	  //cout << "Debug: no "<<Col->getTypeName()<<" Collection is found!" << endl;
          return false;
      }
    return true;

}

/*
        string initString = "M:3,S-1:3,I:9,J:9,K-1:6";          //Need to verify
        isohitcoll->parameters().setValue(LCIO::CellIDEncoding, initString);

        LCFlagImpl flag;
        flag.setBit(LCIO::CHBIT_LONG);
        flag.setBit(LCIO::CHBIT_ID1);
        flag.setBit(LCIO::RCHBIT_ENERGY_ERROR);
        isohitcoll->setFlag(flag.getFlag());
*/


bool SLCIORdr::getHits(LCCollection* SimCol, CellIDDecoder<SimCalorimeterHit>& idDecoderSim, LCCollection* DigiCol, CellIDDecoder<CalorimeterHit>& idDecoderDigi,float x, float y, float z, std::vector<double>& Hits, float& cell_x, float& cell_y, float& cell_z){
      if(SimCol && DigiCol)
      { 
         float min_dist = 1e7;
         int first_index = -1;
         int StaveNum = 0 ;
         int MNum     = 0 ;
         int INum     = 0 ;
         int ZNum     = 0 ;
         int LayerNum = 0 ;
         int aStaveNum = 0 ;
         int aMNum     = 0 ;
         int aINum     = 0 ;
         int aZNum     = 0 ;
         int aLayerNum = 0 ;
         int NHitSim = SimCol->getNumberOfElements();
         int NHitDigi= DigiCol->getNumberOfElements();
         float Hits_E[29][31][31] = {0}; 
         for( int i=0; i<NHitSim; i++ )
         {
             EVENT::SimCalorimeterHit* in = dynamic_cast<EVENT::SimCalorimeterHit*>(SimCol->getElementAt(i));
             //IMPL::CalorimeterHitImpl* in = dynamic_cast<IMPL::CalorimeterHitImpl*>(Col->getElementAt(i));
             LayerNum = idDecoderSim(in)["K-1"];
             if(LayerNum!=0) continue;
             float tmp_dist =  (in->getPosition()[0]-x)*(in->getPosition()[0]-x)+(in->getPosition()[1]-y)*(in->getPosition()[1]-y)+(in->getPosition()[2]-z)*(in->getPosition()[2]-z); 
             if( tmp_dist < min_dist )
             {
                 min_dist = tmp_dist;
                 first_index = i ;
             }
         }
         if (first_index!=-1)
         {
             //IMPL::CalorimeterHitImpl* in = dynamic_cast<IMPL::CalorimeterHitImpl*>(Col->getElementAt(first_index));
             EVENT::SimCalorimeterHit* in = dynamic_cast<EVENT::SimCalorimeterHit*>(SimCol->getElementAt(first_index));
             StaveNum = idDecoderSim(in)["S-1"];// from 0 - 7
             MNum     = idDecoderSim(in)["M"  ];// from 1 - 5
             ZNum     = idDecoderSim(in)["J"  ];// from 0 - 89
             INum     = idDecoderSim(in)["I"  ];// from 0 - 169 for layer 0, then decrease for higher layer
             LayerNum = idDecoderSim(in)["K-1"];// from 0 - 28
             cell_x = in->getPosition()[0] ; 
             cell_y = in->getPosition()[1] ; 
             cell_z = in->getPosition()[2] ; 


             for( int i=0; i<NHitDigi; i++ )
             {
                 IMPL::CalorimeterHitImpl* a_hit = dynamic_cast<IMPL::CalorimeterHitImpl*>(DigiCol->getElementAt(i));
                 aStaveNum = idDecoderDigi(a_hit)["S-1"];// from 0 - 7
                 aMNum     = idDecoderDigi(a_hit)["M"  ];// from 1 - 5
                 aZNum     = idDecoderDigi(a_hit)["J"  ];// from 0 - 89
                 aINum     = idDecoderDigi(a_hit)["I"  ];// from 0 - 169 for layer 0, then decrease for higher layer
                 aLayerNum = idDecoderDigi(a_hit)["K-1"];// from 0 - 28
                 //if((aINum==74 || aINum==75 || aINum==73) && (aLayerNum<=6) )cout<<"Hit "<<i <<",CellID0="<<a_hit->getCellID0()<<",StaveNum="<<aStaveNum<<",MNum="<<aMNum<<",INum="<<aINum<<",ZNum="<<aZNum<<",LayerNum="<<aLayerNum<<", x="<<a_hit->getPosition()[0]<<",y="<<a_hit->getPosition()[1]<<",z="<<a_hit->getPosition()[2]<<"\n";
                 if(StaveNum != aStaveNum) continue; // for same Stave only
                 int tmp_L  = aLayerNum; 
                 int tmp_dZ = aZNum + aMNum*90 - (MNum*90 + ZNum) + 15;//should be 0 - 30 
                 int shift = tmp_L != 0 ? int((tmp_L-0.1)/2.0)+1 : 0 ; // due to the I is shifted for each layer: shift 1 for L 1 and 2, shift 2 for L 3 and 4, shift 3 for L 5 and 6 and so on 
                 int tmp_dI = aINum - INum + shift + 15; //should be 0 - 30
                 //std::cout<<"tmp_L="<< aLayerNum<<",tmp_dZ="<<tmp_dZ<<",tmp_dI="<<tmp_dI<<std::endl;
                 if(tmp_dZ >=0 && tmp_dZ <=30 && tmp_dI >=0 && tmp_dI <=30) {Hits_E[tmp_L][tmp_dZ][tmp_dI] = a_hit->getEnergy(); } 
             }
             
             for (int L=0; L<29; L++)
             {
                 for (int dZ=0; dZ<31; dZ++)
                 {
                     for (int dI=0; dI<31; dI++)
                     {
                         //std::cout<<"Hits_E="<< Hits_E[L][dZ][dI] <<std::endl;
                         Hits.push_back(Hits_E[L][dZ][dI]);
                     }
                  
                 }
             
             } 
             
         
         }

       //   cout<<"ymax="<<ymax<<",ymin="<<ymin<<"Lmax="<<Lmax <<",Lmin="<<Lmin<<",Mmax="<<Mmax<<",Mmin="<<Mmin<<",Imax="<<Imax<<",Imin="<<Imin<<",Jmax="<<Jmax<<",Jmin="<<Jmin<<",Smax="<<Smax<<",Smin="<<Smin<<"\n";
      }
      else {
          std::cout<<"can't find col"<<std::endl;
          return false;
      }
    return true;

}

bool SLCIORdr::getHcalHits(LCCollection* DigiCol, const float& x, const float& y, const float& z, std::vector<double>& Hits, const float& r_min, const float& r_max, const int& r_bin, const float& z_half, const int& z_bin, const float& phi_half, const int& phi_bin){
      //float tmp_r_min = 1e5 ;
      //float tmp_r_max = 0   ;
      float r_interval = (r_max - r_min)/ r_bin ;
      float z_interval = z_half/ z_bin ;
      float phi_interval = phi_half/ phi_bin ;
      const int z_bins = 2*z_bin;
      const int phi_bins = 2*phi_bin;
      if(DigiCol){
         int NHitDigi= DigiCol->getNumberOfElements();
         float Hits_E[r_bin][z_bins][phi_bins] = {0}; 
         for( int i=0; i<NHitDigi; i++ )
         {
             IMPL::CalorimeterHitImpl* hit = dynamic_cast<IMPL::CalorimeterHitImpl*>(DigiCol->getElementAt(i));
             float hit_r = sqrt( hit->getPosition()[0]*hit->getPosition()[0] + hit->getPosition()[1]*hit->getPosition()[1] );
             //if(hit_r > tmp_r_max) tmp_r_max = hit_r;
             //if(hit_r < tmp_r_min) tmp_r_min = hit_r;
             int tmp_r_bin = int ( ( hit_r - r_min ) / r_interval ) ; 
             if(tmp_r_bin<0 || tmp_r_bin>=r_bin) continue;
             int tmp_z_bin = int( (hit->getPosition()[2] - z) / z_interval ) + z_bin ;
             if(tmp_z_bin<0 || tmp_z_bin>=(2*z_bin)) continue;
             float real_phi = getPhi(x, y); 
             float hit_phi  = getPhi(hit->getPosition()[0], hit->getPosition()[1]); 
             int tmp_phi_bin = int( (hit_phi - real_phi) / phi_interval ) + phi_bin ;
             if(tmp_phi_bin<0 || tmp_phi_bin>=(2*phi_bin)) continue;
             Hits_E[tmp_r_bin][tmp_z_bin][tmp_phi_bin] = Hits_E[tmp_r_bin][tmp_z_bin][tmp_phi_bin] + hit->getEnergy(); 
         }
         //std::cout<<"tmp_r_max="<<tmp_r_max<<",tmp_r_min="<<tmp_r_min<<std::endl;
         for (int L=0; L<r_bin; L++)
         {
             for (int dZ=0; dZ<(2*z_bin); dZ++)
             {
                 for (int dI=0; dI<(2*phi_bin); dI++)
                 {
                     Hits.push_back(Hits_E[L][dZ][dI]);
                 }
             }
         }
      }

      else {
          return false;
      }
    return true;
}

bool SLCIORdr::clear(){
    m_HitCm_x = 0;
    m_HitCm_y = 0;
    m_HitCm_z = 0;
    vector<double>().swap(m_HitFirst_x      );
    vector<double>().swap(m_HitFirst_y      );
    vector<double>().swap(m_HitFirst_z      );
    vector<double>().swap(m_phi_rotated     );
    vector<double>().swap(m_HitFirst_vtheta );
    vector<double>().swap(m_HitFirst_vphi   );

    m_HitEn_tot = 0;
    vector<double>().swap(m_mc_vertexX);
    vector<double>().swap(m_mc_vertexY);
    vector<double>().swap(m_mc_vertexZ);
    vector<double>().swap(m_mc_Px     );
    vector<double>().swap(m_mc_Py     );
    vector<double>().swap(m_mc_Pz     );
    vector<double>().swap(m_mc_M      );
    vector<double>().swap(m_Hit_x     );
    vector<double>().swap(m_Hit_y     );
    vector<double>().swap(m_Hit_z     );
    vector<double>().swap(m_Hit_E     ); 
    vector<double>().swap(m_HcalHit_x     );
    vector<double>().swap(m_HcalHit_y     );
    vector<double>().swap(m_HcalHit_z     );
    vector<double>().swap(m_HcalHit_E     ); 
    vector<double>().swap(m_mc_Charge  );
    vector<int>  ().swap(m_mc_Pdg      );
    vector<int>  ().swap(m_mc_genStatus); 
    vector<int>  ().swap(m_mc_simStatus); 
    vector<int>  ().swap(m_mc_np       ); 
    vector<int>  ().swap(m_mc_nd       ); 
    vector<int>  ().swap(m_Hit_ID0     );
    vector<int>  ().swap(m_Hit_ID1     );
    vector<int>  ().swap(m_ID_S     );
    vector<int>  ().swap(m_ID_M     );
    vector<int>  ().swap(m_ID_I     );
    vector<int>  ().swap(m_ID_J     );
    vector<int>  ().swap(m_ID_K     );
    vector<double>().swap(m_mc_pHitx       );
    vector<double>().swap(m_mc_pHity       );
    vector<double>().swap(m_mc_pHitz       );
    vector<double>().swap(m_mc_pHit_theta  );
    vector<double>().swap(m_mc_pHit_phi    );
    vector<double>().swap(m_mc_pHit_rotated);
    vector<double>().swap(m_mc_pHit_dz);
    vector<double>().swap(m_mc_pHit_dy);
    vector< vector<double> >().swap(m_Hits);
    vector< vector<double> >().swap(m_HcalHits);


    return true;
}

bool SLCIORdr::isEnd(){
return false;
}

bool SLCIORdr::configure(){
return true;
}

bool SLCIORdr::finish(){
file_out->cd();
tree_out->Write();
file_out->Close();
return true;
}
