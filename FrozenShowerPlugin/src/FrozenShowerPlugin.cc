
// ********************************************************
// *            Mokka    http://mokka.in2p3.fr            *
// *    -- A Detailed Geant4 Simulation for the ILC --    *
// ********************************************************
//
// $Id: FrozenShowerPlugin.cc,v 1.3 2008/02/08 14:01:24 adrian Exp $
// $Name: mokka-06-06 $
//
#include <string> 
#include <vector>
#include "FrozenShowerPlugin.hh"
#include "PluginManager.hh"
#include "UserInit.hh"
#include "Control.hh"
#include "CGADefs.h"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "CalHit.hh"
#include "VSensitiveDetector.hh"
#include "SEcalSD02.hh"
#include "Encoder32.hh"
#include <math.h>
#include <sstream>
#include <G4Navigator.hh>
#include <ctime>
#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
using namespace AIDA;
#endif
#include "CLHEP/Random/RandGauss.h"
//#define FS_DEBUG  
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
// bins in x, theta,phi,mom
// save found sp for energy shift check
// using Z, X, Y SYMMETRY FOR Ecal Endcap part
// add Ecal Endcap part
// add Energy tunning for barrel
// hit position from cell center
//
template <class Type>
Type stringToNum(const string& str)
{
        istringstream iss(str);
        Type num;
        iss >> num;
        return num;
}


namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}


INITPLUGIN(FrozenShowerPlugin, "FrozenShowerPlugin")

void FrozenShowerPlugin::Init()
{
  std::cout << "#######################################################################################################" << std::endl;
  std::cout << "Using FrozenShowerPlugin " << std::endl 
	    << "  Particles reaching the conditio will be stopped and not tracked any further." 
	    << std::endl ;
  //### ##############
  //m_engine = new HepRandomEngine();
  //m_distribution = new RandGauss(m_engine, 0, 0.1);
  //### ##############
  doSaveSP    = false ;//check the phi cut part!!
  checkTime   = true;
  saveNotFind = false;
  saveFind    = false;
  doEnergyTunning = false;
  apply_Region = 0; // 0 for barrel, 1 for endcap, 2 for barrel + endcap

  theEncoder= new Encoder32();
  width_M = 2350*2/5;
  width_J = width_M/90;

  PI = acos(-1);
  m_x_min = 1850;
  m_x_max = 2000;//mm
  abs_endcap_z_min = 2450;//mm
  abs_endcap_z_max = 2630;//mm
  abs_endcap_inner_x_min = 600;//mm
  abs_endcap_inner_y_min = 600;//mm
  abs_endcap_outer_r_max = 1800;//mm
  Endcap_I_index_max = 161 ;
  Endcap_I_index_min = 0 ;
  Endcap_J_index_max = 241 ;
  Endcap_J_index_min = 0 ;
  Endcap_x_max = 2100;//mm
  Endcap_y_max = 2100;//mm
  Endcap_I_index_bin = (Endcap_x_max - 400.0)/Endcap_I_index_max ;
  Endcap_J_index_bin = (Endcap_y_max + 400.0)/Endcap_J_index_max ;
  std::cout<<"Endcap_I_index_bin="<<Endcap_I_index_bin<<",Endcap_J_index_bin="<<Endcap_J_index_bin<<std::endl;
  m_theta_min = 50;
  m_theta_max = 90;
  //m_phi_min = -10;
  //m_phi_max = 10;
  m_phi_min = -20;
  m_phi_max = 20;
  //m_Mom_min = 100;//MeV
  m_Mom_min = 150;//MeV
  m_Mom_max = 1000;//MeV
  //m_Mom_max = 2000;//MeV
  m_res = 0.15 ;
  m_global_scale = 1.003; //final global E scale value
  //m_Mom_max = 200;//MeV. check with 500MeV
  //m_Mom_max = 150;//MeV. check with 500MeV

  y_boundry = 500;//mm 765
  z_boundry = 2200;//mm 2350
  
  outCut = 0.02;//2%

  //std::istringstream string_rad((*Control::globalModelParameters)["tracker_region_rmax"]);
  //std::istringstream string_z((*Control::globalModelParameters)["tracker_region_zmax"]);
  //string_rad >> tracking_radius_max;
  //string_z >> tracking_z_max;
  if(doSaveSP) std::cout << "Using save start point model" << std::endl;
  else{

  if (doEnergyTunning){
      TFile f0("/junofs/users/wxfang/CEPC/cepcsoft/0.1.0/Simulation/MokkaC/0.1.1/source/Plugin/FrozenShowerPlugin/src/E_tunning.root","read");
      TGraphErrors* gr0 = (TGraphErrors*)f0.Get("ratio_gr_Phibin_E_mean_v1");
      m_gr_phibin_Emean = new TGraphErrors( *gr0 );
      f0.Close();
      for(unsigned i=0; i< m_gr_phibin_Emean->GetN(); i++)
      {
          std::cout<<"i="<<i<<",x="<<m_gr_phibin_Emean->GetX()[i]<<",y="<<m_gr_phibin_Emean->GetY()[i]<<std::endl;
      }
  }
  /*
  std::string barrel_mc_info_em = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_addNoFind//lib_e-.txt";
  std::string barrel_root_em    = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_addNoFind//lib_e-.root";
  std::string barrel_mc_info_ep = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_addNoFind//lib_e+.txt";
  std::string barrel_root_ep    = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_addNoFind//lib_e+.root";
  std::string endcap_mc_info_em = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_Endcap/lib_e-.txt";
  std::string endcap_root_em    = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_Endcap/lib_e-.root";
  std::string endcap_mc_info_ep = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_Endcap/lib_e+.txt";
  std::string endcap_root_ep    = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_Endcap/lib_e+.root";
  */
  std::string barrel_mc_info_em = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/lib_e-.txt";
  std::string barrel_root_em    = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/lib_e-.root";
  std::string barrel_mc_info_ep = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/lib_e+.txt";
  std::string barrel_root_ep    = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/lib_e+.root";
  /*
  std::string barrel_mc_info_em = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/150_2000MeV/lib_e-.txt";
  std::string barrel_root_em    = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/150_2000MeV/lib_e-.root";
  std::string barrel_mc_info_ep = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/150_2000MeV/lib_e+.txt";
  std::string barrel_root_ep    = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Barrel/150_2000MeV/lib_e+.root";
  */
  std::string endcap_mc_info_em = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Endcap/lib_e-.txt";
  std::string endcap_root_em    = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Endcap/lib_e-.root";
  std::string endcap_mc_info_ep = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Endcap/lib_e+.txt";
  std::string endcap_root_ep    = "/junofs/users/wxfang/CEPC/CEPCOFF/FS_study/fs_library/Endcap/lib_e+.root";
  m_h_mom = new TH1D("h_mom","",int(m_Mom_max-m_Mom_min), m_Mom_min, m_Mom_max);
  m_h_res = new TH1D("h_res","",int(m_Mom_max-m_Mom_min), m_Mom_min, m_Mom_max);
  if(apply_Region==0){
      long n_tot_em = 0;
      long n_tot_ep = 0;
      processTxt( barrel_mc_info_em, m_em_x_barrel, m_em_px_barrel, m_em_py_barrel, m_em_pz_barrel, true, m_em_barrel_x_theta_phi_mom_vec, n_tot_em);
      processTxt( barrel_mc_info_ep, m_ep_x_barrel, m_ep_px_barrel, m_ep_py_barrel, m_ep_pz_barrel, true, m_ep_barrel_x_theta_phi_mom_vec, n_tot_ep);
      createMomHist(m_h_mom, m_em_px_barrel, m_em_py_barrel, m_em_pz_barrel, n_tot_em);
      createMomHist(m_h_mom, m_ep_px_barrel, m_ep_py_barrel, m_ep_pz_barrel, n_tot_ep);
      printHist(m_h_mom);
      createResHist(m_h_mom, m_h_res, m_res);
      printHist(m_h_res);
      processMapping(m_barrel_x_chain_map, "e-", barrel_root_em);
      processMapping(m_barrel_x_chain_map, "e+", barrel_root_ep);
  }
  /*
  else if(apply_Region==1){
      processTxt( endcap_mc_info_em, m_em_z_endcap, m_em_px_endcap, m_em_py_endcap, m_em_pz_endcap, false, m_em_z_endcap_mom_vec);
      processTxt( endcap_mc_info_ep, m_ep_z_endcap, m_ep_px_endcap, m_ep_py_endcap, m_ep_pz_endcap, false, m_ep_z_endcap_mom_vec);
      processMapping(m_endcap_z_chain_map, "e-", endcap_root_em);
      processMapping(m_endcap_z_chain_map, "e+", endcap_root_ep);
  }
  else if(apply_Region==2){
      processTxt( barrel_mc_info_em, m_em_x_barrel, m_em_px_barrel, m_em_py_barrel, m_em_pz_barrel, true, m_em_barrel_x_theta_phi_mom_vec);
      processTxt( barrel_mc_info_ep, m_ep_x_barrel, m_ep_px_barrel, m_ep_py_barrel, m_ep_pz_barrel, true, m_ep_barrel_x_theta_phi_mom_vec);
      processTxt( endcap_mc_info_em, m_em_z_endcap, m_em_px_endcap, m_em_py_endcap, m_em_pz_endcap, false, m_em_z_endcap_mom_vec);
      processTxt( endcap_mc_info_ep, m_ep_z_endcap, m_ep_px_endcap, m_ep_py_endcap, m_ep_pz_endcap, false, m_ep_z_endcap_mom_vec);
      processMapping(m_barrel_x_chain_map, "e-", barrel_root_em);
      processMapping(m_barrel_x_chain_map, "e+", barrel_root_ep);
      processMapping(m_endcap_z_chain_map, "e-", endcap_root_em);
      processMapping(m_endcap_z_chain_map, "e+", endcap_root_ep);
  }
  */
  else{
      std::cout<<"Error, wrong Region="<<apply_Region<<",stop!"<<std::endl;
      throw;
  }
  }//else
  if(checkTime && apply_Region==0)// only for barrel now
  {
      std::cout << "Check time..."<<std::endl;
      std::vector<float>* hits_x;
      std::vector<float>* hits_y;
      std::vector<float>* hits_z;
      std::vector<float>* hits_E;
      int  tag;
      std::vector<int>* hits_I;
      std::vector<int>* hits_K;
      for(unsigned int i = 0; i< 1000000; i++)
      {
          m_barrel_x_chain_map[11]->GetEntry(i);
          hits_x = m_ECAL_Hit_x;
          hits_y = m_ECAL_Hit_y;
          hits_z = m_ECAL_Hit_z;
          hits_E = m_ECAL_Hit_E;
          hits_I = m_ECAL_ID_I;
          hits_K = m_ECAL_ID_K;
          tag    = m_Tag;
          m_check_em_hits_x.push_back(*hits_x);
          m_check_em_hits_y.push_back(*hits_y);
          m_check_em_hits_z.push_back(*hits_z);
          m_check_em_hits_E.push_back(*hits_E);
          m_check_em_hits_I.push_back(*hits_I);
          m_check_em_hits_K.push_back(*hits_K);
          m_check_em_tag   .push_back(tag    );
          m_barrel_x_chain_map[-11]->GetEntry(i);
          hits_x = m_ECAL_Hit_x;
          hits_y = m_ECAL_Hit_y;
          hits_z = m_ECAL_Hit_z;
          hits_E = m_ECAL_Hit_E;
          hits_I = m_ECAL_ID_I;
          hits_K = m_ECAL_ID_K;
          tag    = m_Tag;
          m_check_ep_hits_x.push_back(*hits_x);
          m_check_ep_hits_y.push_back(*hits_y);
          m_check_ep_hits_z.push_back(*hits_z);
          m_check_ep_hits_E.push_back(*hits_E);
          m_check_ep_hits_I.push_back(*hits_I);
          m_check_ep_hits_K.push_back(*hits_K);
          m_check_ep_tag   .push_back(tag    );
      }
  }

  //int Depth[29] = {1850, 1857, 1860, 1868, 1871, 1878, 1881, 1889, 1892, 1899, 1902, 1910, 1913, 1920, 1923, 1931, 1934, 1941, 1944, 1952, 1957, 1967, 1972, 1981, 1986, 1996, 2001, 2011, 2016};
  //for(unsigned int i=0; i<29; i++) m_layer.push_back(Depth[i]);

  //read_SMIJK_ID_x_y_z("/junofs/users/wxfang/CEPC/CEPCOFF/ApplyGan/src/apply/SMIJK_ID_x_y_z_new.txt", SMIJK_map_ID_x_y_z);

  G4String out_name =  Control::GetControl()->LCIOFILENAME ;
  out_name.replace(out_name.find(".slcio"), sizeof(".slcio") - 1, ".root");
  std::cout << " OutName= "<<out_name << std::endl;
  //file_out = new TFile("/junofs/users/wxfang/CEPC/cepcsoft/0.1.0/Simulation/MokkaC/0.1.1/source/Plugin/FrozenShowerPlugin/src/check.root","RECREATE");
  file_out = new TFile(out_name.c_str(),"RECREATE");
  tree_out = new TTree("evt","tree");
  tree_out->Branch("m_return" , &m_return      );
  tree_out->Branch("m_s_pid"  , &m_s_pid       );
  tree_out->Branch("m_s_x"    , &m_s_x         );
  tree_out->Branch("m_s_y"    , &m_s_y         );
  tree_out->Branch("m_s_z"    , &m_s_z         );
  tree_out->Branch("m_s_px"   , &m_s_px        );
  tree_out->Branch("m_s_py"   , &m_s_py        );
  tree_out->Branch("m_s_pz"   , &m_s_pz        );
  tree_out->Branch("m_point_x" , &m_point_x   );
  tree_out->Branch("m_point_y" , &m_point_y   );
  tree_out->Branch("m_point_z" , &m_point_z   );
  tree_out->Branch("m_mom_x"   , &m_mom_x   );
  tree_out->Branch("m_mom_y"   , &m_mom_y   );
  tree_out->Branch("m_mom_z"   , &m_mom_z   );
  tree_out->Branch("m_pid"     , &m_pid     );
  tree_out->Branch("m_pass0"   , &m_pass0   );
  tree_out->Branch("m_pass1"   , &m_pass1   );
  tree_out->Branch("m_pass2"   , &m_pass2   );
  tree_out->Branch("m_noFind_point_x" , &m_noFind_point_x   );
  tree_out->Branch("m_noFind_point_y" , &m_noFind_point_y   );
  tree_out->Branch("m_noFind_point_z" , &m_noFind_point_z   );
  tree_out->Branch("m_noFind_mom_x"   , &m_noFind_mom_x   );
  tree_out->Branch("m_noFind_mom_y"   , &m_noFind_mom_y   );
  tree_out->Branch("m_noFind_mom_z"   , &m_noFind_mom_z   );
  tree_out->Branch("m_noFind_pid"     , &m_noFind_pid     );

  tree_out->Branch("m_Find_point_x" , &m_Find_point_x   );
  tree_out->Branch("m_Find_point_y" , &m_Find_point_y   );
  tree_out->Branch("m_Find_point_z" , &m_Find_point_z   );
  tree_out->Branch("m_Find_mom_x"   , &m_Find_mom_x   );
  tree_out->Branch("m_Find_mom_y"   , &m_Find_mom_y   );
  tree_out->Branch("m_Find_mom_z"   , &m_Find_mom_z   );
  tree_out->Branch("m_Find_sum_hit_E", &m_Find_sum_hit_E   );
  tree_out->Branch("m_Find_sum_hit_e", &m_Find_sum_hit_e   );
  tree_out->Branch("m_Find_pid"     , &m_Find_pid     );
  tree_out->Branch("m_FindLib_x"       , &m_FindLib_x   );
  tree_out->Branch("m_FindLib_mom_x"   , &m_FindLib_mom_x   );
  tree_out->Branch("m_FindLib_mom_y"   , &m_FindLib_mom_y   );
  tree_out->Branch("m_FindLib_mom_z"   , &m_FindLib_mom_z   );

}

void FrozenShowerPlugin::Exit()
{
  std::cout << " FrozenShowerPlugin::Exit " << std::endl;
  file_out->cd();
  tree_out->Write();
  file_out->Close();
}

void FrozenShowerPlugin::BeginOfRunAction(const G4Run *run)
{
  //  std::cout << " in FrozenShowerPlugin::BeginOfRunAction " << run->GetRunID() << std::endl;
  //run = 0; // avoid a compiler warning
  int a=0;
}

void FrozenShowerPlugin::EndOfRunAction(const G4Run *run)
{
  //  std::cout << " in FrozenShowerPlugin::EndOfRunAction " << run->GetRunID() << std::endl;
  //run = 0; // avoid a compiler warning
  int a=0;
}

void FrozenShowerPlugin::BeginOfEventAction(const G4Event *evt)
{
  //  std::cout << " in FrozenShowerPlugin::BeginOfEventAction " << evt->GetEventID() << std::endl;
  //evt = 0; // avoid a compiler warning
  m_edep_sum = 0;
  m_pass0 = 0;
  m_pass1 = 0;
  m_pass2 = 0;
  int a=0;
  vector<int>().swap(m_return);
  vector<int>().swap(m_s_pid);
  vector<double>().swap(m_s_x );
  vector<double>().swap(m_s_y );
  vector<double>().swap(m_s_z );
  vector<double>().swap(m_s_px);
  vector<double>().swap(m_s_py);
  vector<double>().swap(m_s_pz);

  vector<double>().swap(m_point_x);
  vector<double>().swap(m_point_y);
  vector<double>().swap(m_point_z);
  vector<double>().swap(m_mom_x  );
  vector<double>().swap(m_mom_y  );
  vector<double>().swap(m_mom_z  );
  vector<int>()   .swap(m_pid  );
  vector<double>().swap(m_noFind_point_x);
  vector<double>().swap(m_noFind_point_y);
  vector<double>().swap(m_noFind_point_z);
  vector<double>().swap(m_noFind_mom_x  );
  vector<double>().swap(m_noFind_mom_y  );
  vector<double>().swap(m_noFind_mom_z  );
  vector<int>()   .swap(m_noFind_pid  );

  vector<double>().swap(m_Find_point_x);
  vector<double>().swap(m_Find_point_y);
  vector<double>().swap(m_Find_point_z);
  vector<double>().swap(m_Find_mom_x  );
  vector<double>().swap(m_Find_mom_y  );
  vector<double>().swap(m_Find_mom_z  );
  vector<double>().swap(m_Find_sum_hit_E);
  vector<double>().swap(m_Find_sum_hit_e);
  vector<int>()   .swap(m_Find_pid  );
  vector<double>().swap(m_FindLib_x);
  vector<double>().swap(m_FindLib_mom_x);
  vector<double>().swap(m_FindLib_mom_y);
  vector<double>().swap(m_FindLib_mom_z);


  //########################################
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  theSD_barrel = dynamic_cast<SEcalSD02*> ( SDman->FindSensitiveDetector("EcalBarrelSilicon") ); 
  if(!theSD_barrel) {std::cout<<"NotFoundSD EcalBarrelSilicon"<<std::endl; throw ;}
  m_barrel_SD_mag = theSD_barrel->CellDim.mag();
  barrel_CalCollection = theSD_barrel->NormalCalCollection;

}

void FrozenShowerPlugin::EndOfEventAction(const G4Event *evt)
{
  //std::cout << " in FrozenShowerPlugin::EndOfEventAction " << evt->GetEventID() << std::endl;
  //std::cout << " sum of edep=" << m_edep_sum << std::endl;
  //evt = 0; // avoid a compiler warning
  int a=0;
  if(checkTime==false) tree_out->Fill();

}

void FrozenShowerPlugin::PreUserTrackingAction(const G4Track *trk)
{
  //  std::cout << " in FrozenShowerPlugin::PreUserTrackingAction " << trk->GetTrackID() << std::endl;

  //trk = 0; // avoid a compiler warning
  int a=0;
}

void FrozenShowerPlugin::PostUserTrackingAction(const G4Track *trk)
{
  //  std::cout << " in FrozenShowerPlugin::PostUserTrackingAction " << trk->GetTrackID() << std::endl;
  //trk = 0; // avoid a compiler warning
  int a=0;
}

void FrozenShowerPlugin::UserSteppingAction(const G4Step *step)
{
//G4cout << "FS, sp x="<<step->GetPreStepPoint()->GetPosition()[0]<<", y="<<step->GetPreStepPoint()->GetPosition()[1]<<", z="<<step->GetPreStepPoint()->GetPosition()[2]<<",id="<<step->GetTrack()->GetTrackID() << G4endl;
int pid = 0;
double x = 0;
double y = 0;
double z = 0;
double px = 0;
double py = 0;
double pz = 0;
int a = makeHitCollection(step, pid, x, y, z, px, py, pz);
m_return.push_back(a);
//if(a==12 || a==13) G4cout << "FS, kill track"<<G4endl;
//else G4cout << "FS, not effect"<<G4endl;
/*
m_s_pid .push_back(pid);
m_s_x   .push_back(x);
m_s_y   .push_back(y);
m_s_z   .push_back(z);
m_s_px  .push_back(px);
m_s_py  .push_back(py);
m_s_pz  .push_back(pz);
*/
}

int FrozenShowerPlugin::makeHitCollection(const G4Step *step, int & _pid, double & _x, double & _y, double & _z, double & _px, double & _py, double & _pz)
{
  //int check = quireIndex(10, 4, "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0514//lib_e-.txt", m_index_em) ;
  //return 0 ;
  //std::cout << " in FrozenShowerPlugin::track " << step->GetTrack()->GetTrackID() << std::endl;
  //t0 = clock();
  G4Track* theTrack  =  step->GetTrack(); 
  const G4ParticleDefinition* def = theTrack->GetDefinition();
  G4int pdgCode=def->GetPDGEncoding();
  const G4ThreeVector &postMom = step->GetPostStepPoint()->GetMomentum();
  float Mom  = sqrt(postMom[0]*postMom[0] + postMom[1]*postMom[1] + postMom[2]*postMom[2]);
  //if( abs(pdgCode) == 22 and Mom < 1 ) {theTrack->SetTrackStatus(fStopAndKill); return 0;}
  if( abs(pdgCode) != 11 ) return 0;
  if( Mom < m_Mom_min || Mom > m_Mom_max ) return 0;

  const G4ThreeVector &postPos = step->GetPostStepPoint()->GetPosition();
  double radius_sqrd = postPos[0]*postPos[0] + postPos[1]*postPos[1];
  /////////////////////////// 
  /*
  if( radius_sqrd > 1850*1850 && abs(postPos[2]) < 2350  ){
      float phi    = getPhi(postPos[0], postPos[1]);
      int part     = partition(phi);// phi should < 360 and > 0
      if (part == 1){ if ( (18<phi && phi<22.5) || (337.5 < phi && phi < 342) return 2;}
      else{
          float rotated = (part-1)*45;
          float tmp_phi = phi-rotated;
          if(abs(tmp_phi>18)) return 2;//gap
      }
  }
  */
  /////////////////////////// 
  bool pass_barrel = true;
  if( radius_sqrd < m_x_min*m_x_min || abs(postPos[2]) > z_boundry  ) pass_barrel=false;
  bool pass_endcap = true;
  if( sqrt(radius_sqrd) > abs_endcap_outer_r_max || ( abs(postPos[0]) < abs_endcap_inner_x_min && abs(postPos[1]) < abs_endcap_inner_y_min ) || abs(postPos[2]) > abs_endcap_z_max || abs(postPos[2]) < abs_endcap_z_min ) pass_endcap=false;
  if     (apply_Region == 0 ){if(pass_barrel == false) return 2;}
  else if(apply_Region == 1 ){if(pass_endcap == false) return 2;}
  else if(apply_Region == 2 ){if(pass_barrel == false && pass_endcap == false) return 2;}
  else return 2 ;
  ////////////////// test speed //////
  //theTrack->SetTrackStatus(fStopAndKill);
  //return 111;
  ///////////////////////////////////
  int Stave; 
  float new_x, new_y, new_px, new_py, rotated ;
  if(pass_barrel){
      if( int(postPos[2]+2350)%188 < 10 || int(postPos[2]+2350)%188 > 178 ) return 3 ; //check gap
      int part = partition(postPos[0],postPos[1]);
      if( part == -1 ) return 3 ;
      Stave = part>=3 ? (part-2) : (part+6);// correct S is from 1 to 8, S-1 is from 0-7 
      float phi    = getPhi(postPos[0], postPos[1]);
      rotated = (part-1)*45;
      float r = sqrt(radius_sqrd);
      new_x = r*cos((phi-rotated)*PI/180);
      new_y = r*sin((phi-rotated)*PI/180);
      if( new_y > 500 || new_y < -750 || new_x<=m_x_min || new_x>=m_x_max ) return 3; //from not found start points, don't consider edge now  
      //////////////remove some phase space may produce hit in HCAL///////////
      /*
      float tmp_x1 = 1970;// mm
      float tmp_y1 = 100;//MeV
      float tmp_x2 = 1900;// mm
      float tmp_y2 = 1000;//MeV
      float tmp_a = (tmp_y2-tmp_y1)/(tmp_x1*tmp_y2-tmp_x2*tmp_y1);
      float tmp_b = (tmp_x1-tmp_x2)/(tmp_x1*tmp_y2-tmp_x2*tmp_y1);
      if ( (tmp_a*new_x + tmp_b*Mom) > 1 ) return 3;// this effect energy scale but not resolution
      */
      //////////////////////////////////////
      float phi_p  = getPhi(postMom[0], postMom[1]);
      float pt = sqrt(postMom[0]*postMom[0] + postMom[1]*postMom[1]);
      new_px = pt*cos((phi_p-rotated)*PI/180);
      new_py = pt*sin((phi_p-rotated)*PI/180);
      if(new_px < 0 ) return 3; //select  
      if(doSaveSP && part !=1) return 3 ; // only use part 1 (abs(phi) <22.5) for save sp
      //if( abs( acos( new_px/sqrt(new_px*new_px + new_py*new_py) )*180/PI ) > 10 )  return 3 ; //FIXME, BE Careful. check phi
      if( abs( acos( new_px/sqrt(new_px*new_px + new_py*new_py) )*180/PI ) > 20 )  return 3 ; //FIXME, BE Careful. check phi
  }
  else if(pass_endcap){
      if(postPos[2]>0){
          if(postMom[2] < 0) return 3; //select  
      }
      else{
          if(postMom[2] > 0) return 3; //select  
      }
  }
  else return 3;
  /////////////// test speed, step 1//////
  theTrack->SetTrackStatus(fStopAndKill);
  return 111;
  /////////////////////////////////
  //t1 = clock();
  //std::cout << "Hi part1, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  //if(apply_Region == 0 && (int(postPos[2]+2350)%188 < 10 || int(postPos[2]+2350)%188 > 178) ) return -1 ; //check gap
  //if(abs(postPos[2]) < z_boundry ){
  //}
  G4double time = step->GetTrack()->GetGlobalTime() ;
  //int part = partition(postPos[0],postPos[1]);
  //if( part == -1 ) return 2 ;
  //int Stave = part>=3 ? (part-2) : (part+6);// correct S is from 1 to 8, S-1 is from 0-7 
  //float phi    = getPhi(postPos[0], postPos[1]);
  //float rotated = (part-1)*45;
  //float r = sqrt(postPos[0]*postPos[0] + postPos[1]*postPos[1]);
  //float new_x = r*cos((phi-rotated)*PI/180);
  //float new_y = r*sin((phi-rotated)*PI/180);
  //std::cout <<"part="<<part<<"phi="<<phi<<",x="<<postPos[0]<<",y="<<postPos[1]<<  std::endl;
  ///////!!!/////////////
  ///////!!!////////////
  //float cos_theta = postMom[2] / Mom;
  //float theta = acos(cos_theta)*180/PI;//0-180
  //if(new_px < 0 || new_x<1870 || new_x>1997 ) return 3; //select  
  if(doSaveSP)
  {
      //float abs_phi = abs( acos( new_px/sqrt(new_px*new_px + new_py*new_py) )*180/PI );
      //if( pass_barrel && ( abs_phi>20 || abs_phi<10 ) ) return 11 ; //FIXME, BE Careful. check phi
      m_point_x.push_back(postPos[0]);
      m_point_y.push_back(postPos[1]);
      m_point_z.push_back(postPos[2]);
      m_mom_x  .push_back(postMom[0]);
      m_mom_y  .push_back(postMom[1]);
      m_mom_z  .push_back(postMom[2]);
      m_pid    .push_back(pdgCode   );
      theTrack->SetTrackStatus(fStopAndKill);
      return 11;
  }
  long nns = -1; 
  float en_scale = 1 ;
  //t1 = clock();
  //std::cout << "Hi part1, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  if(pdgCode==11)
  {
      float lib_x = 0;
      float lib_px = 0;
      float lib_py = 0;
      float lib_pz = 0;
      //if(pass_barrel) nns = Search_v1(m_x_min         , new_x          , new_px    , new_py    , abs(postMom[2]), m_em_x_barrel, m_em_px_barrel, m_em_py_barrel , m_em_pz_barrel, m_em_barrel_x_theta_phi_mom_vec, en_scale, lib_x, lib_px, lib_py, lib_pz) ;
      //if(pass_barrel) nns = Search_v1(m_x_min , m_theta_min , m_phi_min , m_Mom_min      , new_x          , new_px    , new_py    , abs(postMom[2]), m_em_x_barrel, m_em_px_barrel, m_em_py_barrel , m_em_pz_barrel, m_em_barrel_x_theta_phi_mom_vec, en_scale, lib_x, lib_px, lib_py, lib_pz) ;
      if(pass_barrel) nns = Search_v1(m_x_min , m_theta_min , m_phi_min , m_Mom_min      , new_x          , new_px    , new_py    , abs(postMom[2]), m_em_x_barrel, m_em_px_barrel, m_em_py_barrel , m_em_pz_barrel, m_em_barrel_x_theta_phi_mom_vec, en_scale, lib_x, lib_px, lib_py, lib_pz, m_h_res) ;
      else            nns = Search(abs_endcap_z_min, abs(postPos[2]), abs(postMom[0]), abs(postMom[1]), abs(postMom[2]), m_em_z_endcap, m_em_px_endcap, m_em_py_endcap , m_em_pz_endcap, m_em_z_endcap_mom_vec, en_scale) ;
      if(nns != -1 )
      {
          if(checkTime==false){
                              if(pass_barrel) m_barrel_x_chain_map[11]->GetEntry(nns);
                              else            m_endcap_z_chain_map[11]->GetEntry(nns);
                               //std::cout<<"em nns="<<nns<<std::endl;
                              
                              if(saveFind){
                                  m_Find_point_x.push_back(postPos[0]);
                                  m_Find_point_y.push_back(postPos[1]);
                                  m_Find_point_z.push_back(postPos[2]);
                                  m_Find_mom_x  .push_back(postMom[0]);
                                  m_Find_mom_y  .push_back(postMom[1]);
                                  m_Find_mom_z  .push_back(postMom[2]);
                                  m_Find_pid    .push_back(pdgCode   );
                                  m_FindLib_x       .push_back(lib_x);
                                  m_FindLib_mom_x   .push_back(lib_px);
                                  m_FindLib_mom_y   .push_back(lib_py);
                                  m_FindLib_mom_z   .push_back(lib_pz);
                                  }
                               
                              }
          else
          {
              int c_size = m_check_em_hits_x.size();
              m_ECAL_Hit_x = &(m_check_em_hits_x.at(nns%c_size));   
              m_ECAL_Hit_y = &(m_check_em_hits_y.at(nns%c_size));   
              m_ECAL_Hit_z = &(m_check_em_hits_z.at(nns%c_size));   
              m_ECAL_Hit_E = &(m_check_em_hits_E.at(nns%c_size));   
              m_ECAL_ID_I = &(m_check_em_hits_I .at(nns%c_size));   
              m_ECAL_ID_K = &(m_check_em_hits_K .at(nns%c_size));   
              m_Tag        = m_check_em_tag     .at(nns%c_size) ;   
          }
      }
      else
      {
          if(saveNotFind)
          {
              if(pass_barrel){
              m_noFind_point_x.push_back(new_x);
              m_noFind_point_y.push_back(new_y);
              m_noFind_point_z.push_back(postPos[2]);
              m_noFind_mom_x  .push_back(new_px);
              m_noFind_mom_y  .push_back(new_py);
              m_noFind_mom_z  .push_back(postMom[2]);
              m_noFind_pid    .push_back(pdgCode   );
              }
              else{
              m_noFind_point_x.push_back(postPos[0]);
              m_noFind_point_y.push_back(postPos[1]);
              m_noFind_point_z.push_back(abs(postPos[2]));
              m_noFind_mom_x  .push_back(abs(postMom[0]));
              m_noFind_mom_y  .push_back(abs(postMom[1]));
              m_noFind_mom_z  .push_back(abs(postMom[2]));
              m_noFind_pid    .push_back(pdgCode   );
              }
          }
      return 4;
      }
  }
  else if(pdgCode==-11)
  {
      float lib_x = 0;
      float lib_px = 0;
      float lib_py = 0;
      float lib_pz = 0;
      //if(pass_barrel) nns = Search_v1(m_x_min         , new_x          , new_px    , new_py    , abs(postMom[2]), m_ep_x_barrel, m_ep_px_barrel, m_ep_py_barrel , m_ep_pz_barrel, m_ep_barrel_x_theta_phi_mom_vec, en_scale, lib_x, lib_px, lib_py, lib_pz) ;
      //if(pass_barrel) nns = Search_v1(m_x_min , m_theta_min , m_phi_min , m_Mom_min , new_x          , new_px    , new_py    , abs(postMom[2]), m_ep_x_barrel, m_ep_px_barrel, m_ep_py_barrel , m_ep_pz_barrel, m_ep_barrel_x_theta_phi_mom_vec, en_scale, lib_x, lib_px, lib_py, lib_pz) ;
      if(pass_barrel) nns = Search_v1(m_x_min , m_theta_min , m_phi_min , m_Mom_min , new_x          , new_px    , new_py    , abs(postMom[2]), m_ep_x_barrel, m_ep_px_barrel, m_ep_py_barrel , m_ep_pz_barrel, m_ep_barrel_x_theta_phi_mom_vec, en_scale, lib_x, lib_px, lib_py, lib_pz, m_h_res) ;
      else            nns = Search(abs_endcap_z_min, abs(postPos[2]), abs(postMom[0]), abs(postMom[1]), abs(postMom[2]), m_ep_z_endcap, m_ep_px_endcap, m_ep_py_endcap , m_ep_pz_endcap, m_ep_z_endcap_mom_vec, en_scale) ;
      if(nns != -1 ) 
      {
          if(checkTime==false){if(pass_barrel)m_barrel_x_chain_map[-11]->GetEntry(nns);
                               else           m_endcap_z_chain_map[-11]->GetEntry(nns);
                               //std::cout<<"ep nns="<<nns<<std::endl;
                               if(saveFind){
                                   m_Find_point_x.push_back(postPos[0]);
                                   m_Find_point_y.push_back(postPos[1]);
                                   m_Find_point_z.push_back(postPos[2]);
                                   m_Find_mom_x  .push_back(postMom[0]);
                                   m_Find_mom_y  .push_back(postMom[1]);
                                   m_Find_mom_z  .push_back(postMom[2]);
                                   m_Find_pid    .push_back(pdgCode   );
                                   m_FindLib_x       .push_back(lib_x);
                                   m_FindLib_mom_x   .push_back(lib_px);
                                   m_FindLib_mom_y   .push_back(lib_py);
                                   m_FindLib_mom_z   .push_back(lib_pz);
                                   }
                              }
          else
          {
              int c_size   = m_check_ep_hits_x.size();
              m_ECAL_Hit_x = &(m_check_ep_hits_x.at(nns%c_size));   
              m_ECAL_Hit_y = &(m_check_ep_hits_y.at(nns%c_size));   
              m_ECAL_Hit_z = &(m_check_ep_hits_z.at(nns%c_size));   
              m_ECAL_Hit_E = &(m_check_ep_hits_E.at(nns%c_size));   
              m_ECAL_ID_I = &(m_check_ep_hits_I .at(nns%c_size));   
              m_ECAL_ID_K = &(m_check_ep_hits_K .at(nns%c_size));   
              m_Tag        = m_check_ep_tag     .at(nns%c_size) ;   
          }
      }
      else
      { 
          if(saveNotFind)
          {   
              if(pass_barrel){
              m_noFind_point_x.push_back(new_x);
              m_noFind_point_y.push_back(new_y);
              m_noFind_point_z.push_back(postPos[2]);
              m_noFind_mom_x  .push_back(new_px);
              m_noFind_mom_y  .push_back(new_py);
              m_noFind_mom_z  .push_back(postMom[2]);
              m_noFind_pid    .push_back(pdgCode   );
              }
              else{
              m_noFind_point_x.push_back(postPos[0]);
              m_noFind_point_y.push_back(postPos[1]);
              m_noFind_point_z.push_back(abs(postPos[2]));
              m_noFind_mom_x  .push_back(abs(postMom[0]));
              m_noFind_mom_y  .push_back(abs(postMom[1]));
              m_noFind_mom_z  .push_back(abs(postMom[2]));
              m_noFind_pid    .push_back(pdgCode   );
              }
          }
          return 5;
      }
  }
  //////////////// test speed, step 2//////
  //theTrack->SetTrackStatus(fStopAndKill);
  //return 111;
  /////////////////////////////////
  //t0 = clock();
  //std::cout << "Hi part2, dt="<<(t0-t1) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  //t0 = clock();
  //std::cout << "Hi part2, dt="<<(t0-t1) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  
  //m_pass0 ++ ;
  //if(m_Tag == 2) return 7;// remove when hcal has hit ?
  //m_pass1 ++ ;
  //t1 = clock();
  //std::cout << "Hi 1, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  //t0 = clock();
  //# time consuming part !! ####
  //map<int, TChain*>::iterator it = m_barrel_x_chain_map.find(pdgCode);
  //if( it != m_barrel_x_chain_map.end()){it->second->GetEntry(index);}
  //else {return 8;}
  //std::cout<<"new_x="<<new_x<<",int(new_x)="<<int(new_x)<<",pid="<<_pid<<",x="<<_x<<",y="<<_y<<",z"<<_z<<",_px"<<_px<<",py="<<_py<<",pz="<<_pz<<",bin_theta="<<bin_theta<<",bin_phi="<<bin_phi<<",index="<<index<<",m_Tag->size()="<<m_Tag->size()<<std::endl; 
  //##############
  /*
  t1 = clock();
  std::cout << "Hi map, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  t0 = clock();
  */
  float tot_E = 0;
  float tot_e = 0;
  if(pass_barrel)
  {
     
      G4int n_hit = barrel_CalCollection->entries();
      G4int PDG = pdgCode;
      G4Step* astep = const_cast<G4Step*>(step);
      G4int PID =  Control::GetControl()->GetPIDForCalHit(astep);//for found primary MC particle
      //#####################################
      float en_tunning = doEnergyTunning==true ? getE_scale(postPos[0], postPos[1])*m_global_scale : 1 ;//get e scale from sp
      G4int theSDPiece = 0;
      theSDPiece = ECALBARREL;
      int N_find = 0;
      for(unsigned int i=0; i< m_ECAL_Hit_y->size(); i++)
      {
          float Hit_x     = m_ECAL_Hit_x->at(i) + new_x ;
          int l_index = -1;
          l_index = m_ECAL_ID_K->at(i) + 2;//due to get from K-1, the preshower is 1 then 2, 3,4... for Ecal
          if(l_index==-1) continue;
          float Hit_y     = m_ECAL_Hit_y->at(i) ;
          float Hit_z     = m_ECAL_Hit_z->at(i) ;
          // parity with x-z plane 
          float new_Hit_y = Hit_y ;
          // parity with x-y plane 
          float new_Hit_z = postMom[2] > 0 ? Hit_z : -Hit_z ;
          //moving
          new_Hit_y = new_Hit_y + new_y;//from center to new y
          int I_index = m_ECAL_ID_I->at(i) - int(new_y/10);// new I
          new_Hit_z = new_Hit_z + postPos[2];
          // rotating with z axis 
          float new_Hit_r = sqrt( Hit_x*Hit_x + new_Hit_y*new_Hit_y );
          float real_Hit_z =  new_Hit_z ;
          float new_phi = getPhi(Hit_x, new_Hit_y);
          float real_Hit_x = new_Hit_r*cos( (new_phi+rotated)*PI/180 );
          float real_Hit_y = new_Hit_r*sin( (new_phi+rotated)*PI/180 );
          //std::cout<<"mc_x="<<new_x<<",mc_theta="<<theta<<",mc_mom="<<Mom<<",index="<<index<<",lib_index="<<lib_index << ",real x="<<real_Hit_x<<", real y="<<real_Hit_y<<", real z="<<real_Hit_z<<", en="<<m_ECAL_Hit_E->at(lib_index).at(i)<< std::endl;
          int M_index = int((real_Hit_z+2350)/width_M) + 1 ; // 1 to 5 
          int J_index = int((real_Hit_z+2350 - (M_index-1)*width_M)/width_J) ; // 0 to 89 
          if(I_index <0) I_index = 0;
          else if(I_index > 169) I_index = 169;
          if(M_index < 1) M_index = 1;
          else if(M_index > 5) M_index = 5;
          if(J_index < 0) J_index = 0;
          else if(J_index > 89) J_index = 89;
         
          G4ThreeVector thePosition(real_Hit_x, real_Hit_y, real_Hit_z);   
          G4ThreeVector theCellCenter(0, 0, 0);   
          G4ThreeVector dis(0, 0, 0);   
          int final_S = Stave;
          int final_M = M_index;
          int final_J = J_index;
          int final_I = I_index;
          bool found_id = false;
          bool found_best_id = false;
          float min_dis = 999;
          float secd_min_dis = 999;
      //t0 = clock();
          
          for(int tmp_M = (M_index-1) < 1 ? 1 : M_index-1, max_M = (M_index+1) > 5 ? 5 : M_index+1 ; tmp_M<=max_M; tmp_M++)
          {
              for(int tmp_I = (I_index-10) < 0 ? 0 : I_index-10 , max_I = (I_index+10) > 169 ? 169 : I_index+10 ; tmp_I<=max_I; tmp_I++)
              {
                  for(int tmp_J = 0; tmp_J<=89; tmp_J++)
                  {
                      int tmp_S = Stave;
                      theCellCenter = theSD_barrel->GetCellCenter(theSDPiece, tmp_S, tmp_M, tmp_I, tmp_J, l_index); 
                      dis = theCellCenter - thePosition ;
                      if ( (dis.mag() < min_dis) && (dis.mag() < sqrt(2*m_barrel_SD_mag*m_barrel_SD_mag+ 3*3) ) ) 
                      {
                          secd_min_dis = min_dis;
                          min_dis = dis.mag();
                          //final_S = Stave;
                          final_S = tmp_S;
                          final_M = tmp_M;
                          final_I = tmp_I;
                          final_J = tmp_J;
                          found_id = true;
                          if(min_dis < 5 ) {found_best_id = true; break;}
                      }
                  }
                  if(found_best_id) break;
              }
              if(found_best_id) break;
          }
      //t1 = clock();
      //std::cout << "Hi loop, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
      //t0 = clock();
          if(found_id){
              G4ThreeVector finalCellCenter = theSD_barrel->GetCellCenter(theSDPiece, final_S, final_M, final_I, final_J, l_index);
              cell_ids theIndex             =  theEncoder->encode  (final_S, final_M, final_I, final_J, l_index-1, 0);
              float edep = m_ECAL_Hit_E->at(i)*1000*en_scale*en_tunning;//to MeV
              bool found = false;
              for(G4int i_hit = 0; i_hit< n_hit; i_hit++)
              {
                if((*barrel_CalCollection)[i_hit]->testCell(theIndex)) 
                  {
                  (*barrel_CalCollection)[i_hit]->AddEdep(PID,PDG,edep,time);
                  found = true;
                  break;
                  }
              }
              if(!found)
              {
                  barrel_CalCollection->  insert(new CalHit (theSDPiece,
                			  final_S,
                			  final_M,
                			  final_I,
                			  final_J,
                			  l_index-1,
                			  0,
                			  finalCellCenter(0),
                			  finalCellCenter(1),
                			  finalCellCenter(2),
                			  edep,
                			  PID,
                			  PDG,
                			  time,
                			  theIndex));
              }
          }//found
      }// for hit

  }
  else// for endcap
  {
      //########################################
      G4SDManager* SDman = G4SDManager::GetSDMpointer();
      SEcalSD02* theSD = dynamic_cast<SEcalSD02*> ( SDman->FindSensitiveDetector("EcalEndcapSilicon") ); 
      if(!theSD) {std::cout<<"NotFoundSD EcalEndcapSilicon"<<std::endl; return 11;}
      HitsCollection *CalCollection;
      CalCollection = theSD->NormalCalCollection;
      G4int n_hit = CalCollection->entries();
      G4int PDG = pdgCode;
      G4Step* astep = const_cast<G4Step*>(step);
      G4int PID =  Control::GetControl()->GetPIDForCalHit(astep);//for found primary MC particle
      G4int theSDPiece = 0;
      for(unsigned int i=0; i< m_ECAL_Hit_y->size(); i++)
      {
          int S_index = -1;
          int M_index = -1;
          int I_index = -1;
          int J_index = -1;
          int K_index = -1;
          K_index = m_ECAL_ID_K->at(i) + 2;//due to get from K-1, the preshower is 1 then 2, 3,4... for Ecal
          if(K_index==-1) continue;
          //std::cout<<"K_index="<<K_index<<",Hit_x="<<Hit_x<<std::endl;
          float Hit_x     = postMom[0] > 0 ? m_ECAL_Hit_x->at(i) + postPos[0] : -m_ECAL_Hit_x->at(i) + postPos[0] ;
          float Hit_y     = postMom[1] > 0 ? m_ECAL_Hit_y->at(i) + postPos[1] : -m_ECAL_Hit_y->at(i) + postPos[1] ;
          float Hit_z     = postPos[2] > 0 ? m_ECAL_Hit_z->at(i) + postPos[2] : -m_ECAL_Hit_z->at(i) + postPos[2] ;
          float xy_phi    = getPhi(Hit_x, Hit_y);
          float xy_r      = sqrt(Hit_x*Hit_x + Hit_y*Hit_y);
          if(Hit_z<0){
              theSDPiece =  ECALENDCAPMINUS; 
              float new_Hit_x = 0;
              float new_Hit_y = 0;
              M_index = 0; 
              if     (Hit_x> 400 && Hit_y>-400) { S_index = 0 + 1 ; new_Hit_x = Hit_x; new_Hit_y = Hit_y; } //S-1 
              else if(Hit_x>-400 && Hit_y<-400) { S_index = 1 + 1 ; new_Hit_x = xy_r*cos((xy_phi+90 )*PI/180); new_Hit_y = xy_r*sin((xy_phi+90 )*PI/180); }
              else if(Hit_x<-400 && Hit_y< 400) { S_index = 2 + 1 ; new_Hit_x = xy_r*cos((xy_phi+180)*PI/180); new_Hit_y = xy_r*sin((xy_phi+180)*PI/180); }
              else if(Hit_x< 400 && Hit_y> 400) { S_index = 3 + 1 ; new_Hit_x = xy_r*cos((xy_phi+270)*PI/180); new_Hit_y = xy_r*sin((xy_phi+270)*PI/180); }
              J_index = (new_Hit_y + 400)/Endcap_J_index_bin;
              I_index = (new_Hit_x - 400)/Endcap_I_index_bin;
          }
          else{
              theSDPiece =  ECALENDCAPPLUS; 
              float new_Hit_x = 0;
              float new_Hit_y = 0;
              M_index = 6; 
              if     (Hit_x< 400 && Hit_y<-400) { S_index = 1 + 1 ; new_Hit_x = xy_r*cos((xy_phi+180)*PI/180); new_Hit_y = xy_r*sin((xy_phi+180)*PI/180); }
              else if(Hit_x<-400 && Hit_y>-400) { S_index = 0 + 1 ; new_Hit_x = xy_r*cos((xy_phi+270)*PI/180); new_Hit_y = xy_r*sin((xy_phi+270)*PI/180); }
              else if(Hit_x> 400 && Hit_y< 400) { S_index = 2 + 1 ; new_Hit_x = xy_r*cos((xy_phi+90 )*PI/180); new_Hit_y = xy_r*sin((xy_phi+90 )*PI/180); }
              else if(Hit_x>-400 && Hit_y> 400) { S_index = 3 + 1 ; new_Hit_x = Hit_x; new_Hit_y = Hit_y; }
              J_index = (new_Hit_x + 400)/Endcap_J_index_bin;
              I_index = (new_Hit_y - 400)/Endcap_I_index_bin;
          }
          if(S_index==-1) continue;//out of region 
          G4ThreeVector thePosition(Hit_x, Hit_y, Hit_z);   
          G4ThreeVector theCellCenter(0, 0, 0);   
          G4ThreeVector dis(0, 0, 0);   
          int final_S = S_index;
          int final_M = M_index;
          int final_I = I_index;
          int final_J = J_index;
          int final_K = K_index;
          bool found_id = false;
          bool found_best_id = false;
          float min_dis = 999;
          float secd_min_dis = 999;
      //t0 = clock();
          for(int tmp_I = (I_index-10) < 0 ? 0 : I_index-10 , max_I = (I_index+10) > Endcap_I_index_max ? Endcap_I_index_max : I_index+10 ; tmp_I<max_I; tmp_I++)
          {
              for(int tmp_J = (J_index-10) < 0 ? 0 : J_index-10 , max_J = (J_index+10) > Endcap_J_index_max ? Endcap_J_index_max : J_index+10 ; tmp_J<max_J; tmp_J++)
              {
                  theCellCenter = theSD->GetCellCenter(theSDPiece, S_index, M_index, tmp_I, tmp_J, K_index); 
                  dis = theCellCenter - thePosition ;
                  if ( (dis.mag() < min_dis) && (dis.mag() < theSD->CellDim.mag()*sqrt(2.))) 
                  {
                      secd_min_dis = min_dis;
                      min_dis = dis.mag();
                      final_I = tmp_I;
                      final_J = tmp_J;
                      found_id = true;
                      if(min_dis < 5 ) {found_best_id = true; break;}
                  }
              }
              if(found_best_id) break;
          }
          //std::cout<<"Hit_x="<<Hit_x<<",Hit_y="<<Hit_y<<",Hit_z="<<Hit_z<<",theSDPiece="<<theSDPiece<<",final_S="<<final_S<<",final_M="<<final_M<<",final_I="<<final_I<<",final_J="<<final_J<<",final_K="<<final_K<<std::endl;
          G4ThreeVector finalCellCenter = theSD->GetCellCenter(theSDPiece, final_S, final_M, final_I, final_J, final_K);
          //std::cout<<"finalCellCenter(0)="<<finalCellCenter(0)<<",finalCellCenter(1)="<<finalCellCenter(1)<<",finalCellCenter(2)="<<finalCellCenter(2)<<std::endl;
          if(found_id == false)
          {
              finalCellCenter =  thePosition ;
          }
          cell_ids theIndex =  theEncoder->encode(final_S, final_M, final_I, final_J, final_K-1, 0);// 0 means ?
          //std::cout<<"Endcap, min_dis="<<min_dis<<",secd mis="<<secd_min_dis<<",found_id="<<found_id<<",CellDim.mag()="<<theSD->CellDim.mag()<<",x="<<finalCellCenter(0)<<",y="<<finalCellCenter(1)<<",z="<<finalCellCenter(2)<<",real x="<<Hit_x<<", real y="<<Hit_y<<", real z="<<Hit_z<<",sp x="<<postPos[0]<<",sp y="<<postPos[1]<<",sp z="<<postPos[2]<<",mom x="<<postMom[0]<<",mom y="<<postMom[1]<<",mom z="<<postMom[2]<< std::endl;
      //#####################################
          float edep = m_ECAL_Hit_E->at(i)*1000*en_scale ;//to MeV
          bool found = false;
          for(G4int i_hit = 0; i_hit< n_hit; i_hit++)
          {
            if((*CalCollection)[i_hit]->testCell(theIndex)) 
              {
              (*CalCollection)[i_hit]->AddEdep(PID,PDG,edep,time);
              found = true;
              break;
              }
          }
          if(!found)
          {
              CalCollection->  insert(new CalHit (theSDPiece,
            			  final_S,
            			  final_M,
            			  final_I,
            			  final_J,
            			  final_K-1,
            			  0,
            			  finalCellCenter(0),
            			  finalCellCenter(1),
            			  finalCellCenter(2),
            			  edep,
            			  PID,
            			  PDG,
            			  time,
            			  theIndex));
          }
      }
  }//endcap
  if(saveFind){
      m_Find_sum_hit_E.push_back(tot_E);
      m_Find_sum_hit_e.push_back(tot_e);
  }
  
  theTrack->SetTrackStatus(fStopAndKill);
  //t1 = clock();
  //std::cout << "Hi part3, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t1="<<t1<<",t0="<<t0<<std::endl;
  return 12;      
}

float FrozenShowerPlugin::getE_scale(float x, float y) // from -90 to 90
{
    float phi = getPhi(x, y);
    if(90<phi && phi < 270)     phi = 180 - phi;
    else if(270<phi && phi<360) phi = phi - 360;
    if(phi<m_gr_phibin_Emean->GetX()[0]) return m_gr_phibin_Emean->GetY()[0];
    else if(phi>m_gr_phibin_Emean->GetX()[m_gr_phibin_Emean->GetN()-1]) return m_gr_phibin_Emean->GetY()[m_gr_phibin_Emean->GetN()-1];
    else{
        for(unsigned i=0; i< m_gr_phibin_Emean->GetN()-1; i++)
        {
            if(m_gr_phibin_Emean->GetX()[i]<phi && m_gr_phibin_Emean->GetX()[i+1] > phi)
            {
                return (m_gr_phibin_Emean->GetY()[i] + (phi-m_gr_phibin_Emean->GetX()[i])*(m_gr_phibin_Emean->GetY()[i+1]-m_gr_phibin_Emean->GetY()[i])/(m_gr_phibin_Emean->GetX()[i+1]-m_gr_phibin_Emean->GetX()[i])  );
            }
            //std::cout<<"i="<<i<<",x="<<m_gr_phibin_Emean->GetX()[i]<<",y="<<m_gr_phibin_Emean->GetY()[i]<<std::endl;
        }
    }
    return 1;
}


void FrozenShowerPlugin::getID_x_y_z(const int& S, const int& M, const int& I, const int& J, const int& K, int& id, float& x, float& y, float& z)
{
    vector<string> str_result;
    int result = split(str_result, SMIJK_map_ID_x_y_z[S][M][I][J][K] , "_");
    if(str_result.size()==4)
    {
        id = stringToNum <int> (str_result.at(0));
        x  = stringToNum<float>(str_result.at(1));
        y  = stringToNum<float>(str_result.at(2));
        z  = stringToNum<float>(str_result.at(3));
    }
    else
    {
        id = 0;
        x  = 0;
        y  = 0;
        z  = 0;
    }
}

int FrozenShowerPlugin::split(vector<string>& res, const string& str, const string& delim) {
    if("" == str) return 0;  
    //先将要切割的字符串从string类型转换为char*类型  
    char* strs = new char[str.length() + 1] ;
    strcpy(strs, str.c_str());   
    char* d = new char[delim.length() + 1];  
    strcpy(d, delim.c_str());  
    char* p = strtok(strs, d);  
    while(p) 
    {  
        string s = p; //分割得到的字符串转换为string类型  
        res.push_back(s); //存入结果数组 
        p = strtok(NULL, d);
    }
    delete [] strs;  
    delete [] d;  
    return 1;  
    
}


long FrozenShowerPlugin::Search(const float& low, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector<long> > > & conts, float& en_scale)
{ 
    long index=-1;
    float min_dist = 10;
    float tmp_mom = sqrt(q_px*q_px + q_py*q_py  + q_pz*q_pz);
    int flow_x = 1 ;
    int range_low_x  = int(q_x-low-flow_x) > 0            ? int(q_x-low-flow_x) :            0   ;
    int range_high_x = int(q_x-low+flow_x) < conts.size() ? int(q_x-low+flow_x) : conts.size()-1 ;
    int flow_mom = 5 ;
    int range_low_mom  = int(tmp_mom-m_Mom_min-flow_mom) > 0                  ? int(tmp_mom-m_Mom_min-flow_mom) :                  0   ;
    int range_high_mom = int(tmp_mom-m_Mom_min+flow_mom) < conts.at(0).size() ? int(tmp_mom-m_Mom_min+flow_mom) : conts.at(0).size()-1 ;
    for(int i = range_low_x ; i <= range_high_x; i++ )
    {
        for(int j=range_low_mom ; j < range_high_mom ; j++)
        {
            for(long k=0 ; k < conts.at(i).at(j).size() ; k++)
            {
                long in = conts.at(i).at(j).at(k);
                float pre_cut = sqrt(min_dist);
                if(abs(db_x[in]-q_x) > pre_cut || abs(db_px[in]-q_px) > pre_cut || abs(db_py[in]-q_py) > pre_cut || abs(db_pz[in]-q_pz) > pre_cut ) continue;
                float dist =  (db_x[in]-q_x)*(db_x[in]-q_x) + (db_px[in]-q_px)*(db_px[in]-q_px) + (db_py[in]-q_py)*(db_py[in]-q_py) + (db_pz[in]-q_pz)*(db_pz[in]-q_pz) ;
                if( dist < min_dist ) {min_dist = dist; index = in;}
                //if( min_dist < 1) break;
            }
        }
    }
    en_scale = (index != -1) ? tmp_mom/sqrt(db_px[index]*db_px[index] + db_py[index]*db_py[index] + db_pz[index]*db_pz[index]) : 1 ;
    //if(index!=-1)std::cout<<"L2 vector conts_size="<<conts.size()<<",index="<<index<<",en_scale="<<en_scale<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<",dist="<< min_dist <<std::endl;
    return index; 
}

long FrozenShowerPlugin::Search(const float& low, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector<long> > > & conts, float& en_scale, float& b_x, float& b_px, float& b_py, float& b_pz)
{ 
    long index=-1;
    float min_dist = 10;
    float tmp_mom = sqrt(q_px*q_px + q_py*q_py  + q_pz*q_pz);
    int flow_x = 1 ;
    //int flow_x = 0 ;
    int range_low_x  = int(q_x-low-flow_x) > 0            ? int(q_x-low-flow_x) :            0   ;
    int range_high_x = int(q_x-low+flow_x) < conts.size() ? int(q_x-low+flow_x) : conts.size()-1 ;
    int flow_mom = 5 ;
    int range_low_mom  = int(tmp_mom-m_Mom_min-flow_mom) > 0                  ? int(tmp_mom-m_Mom_min-flow_mom) :                  0   ;
    int range_high_mom = int(tmp_mom-m_Mom_min+flow_mom) < conts.at(0).size() ? int(tmp_mom-m_Mom_min+flow_mom) : conts.at(0).size()-1 ;
    for(int i = range_low_x ; i <= range_high_x; i++ )
    {
        for(int j=range_low_mom ; j < range_high_mom ; j++)
        {
            for(long k=0 ; k < conts.at(i).at(j).size() ; k++)
            {
                long in = conts.at(i).at(j).at(k);
                float pre_cut = sqrt(min_dist);
                if(abs(db_x[in]-q_x) > pre_cut || abs(db_px[in]-q_px) > pre_cut || abs(db_py[in]-q_py) > pre_cut || abs(db_pz[in]-q_pz) > pre_cut ) continue;
                //if(abs(db_px[in]-q_px) > 1 || abs(db_py[in]-q_py) > 1 || abs(db_pz[in]-q_pz) > 1 ) continue; //FIXME
                float dist =  (db_x[in]-q_x)*(db_x[in]-q_x) + (db_px[in]-q_px)*(db_px[in]-q_px) + (db_py[in]-q_py)*(db_py[in]-q_py) + (db_pz[in]-q_pz)*(db_pz[in]-q_pz) ;
                if( dist < min_dist ) {min_dist = dist; index = in;}
                //if( min_dist < 1) break;
            }
        }
    }
    en_scale = (index != -1) ? tmp_mom/sqrt(db_px[index]*db_px[index] + db_py[index]*db_py[index] + db_pz[index]*db_pz[index]) : 1 ;
    //if(index!=-1)std::cout<<"L2 vector conts_size="<<conts.size()<<",index="<<index<<",en_scale="<<en_scale<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<",dist="<< min_dist <<std::endl;
    if(index != -1){b_x=db_x[index]; b_px=db_px[index]; b_py=db_py[index]; b_pz=db_pz[index];}
    return index; 
}

long FrozenShowerPlugin::Search_v1(const float& low, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector<long> > > & conts, float& en_scale, float& b_x, float& b_px, float& b_py, float& b_pz)
{ 
    //float min_dist = 10;
    float tmp_mom = sqrt(q_px*q_px + q_py*q_py  + q_pz*q_pz);

    //float float_value = 0.03*tmp_mom;
    //int flow_x = 2 ;
    //int flow_mom = int(0.05*tmp_mom) ;
    //float float_value = 0.02*tmp_mom;
    //int flow_x = 1 ;
    //int flow_mom = int(0.03*tmp_mom) ;
    float float_value = 0.1*tmp_mom;
    int flow_x = 0 ;
    int flow_mom = int(0.1*tmp_mom) ;
    //float float_value = 0.02*tmp_mom;
    //int flow_x = 0 ;
    //int flow_mom = int(0.03*tmp_mom) ;
    int range_low_x  = int(q_x-low-flow_x) > 0            ? int(q_x-low-flow_x) :            0   ;
    int range_high_x = int(q_x-low+flow_x) < conts.size() ? int(q_x-low+flow_x) : conts.size()-1 ;
    int range_low_mom  = int(tmp_mom-m_Mom_min-flow_mom) > 0                  ? int(tmp_mom-m_Mom_min-flow_mom) :                  0   ;
    int range_high_mom = int(tmp_mom-m_Mom_min+flow_mom) < conts.at(0).size() ? int(tmp_mom-m_Mom_min+flow_mom) : conts.at(0).size()-1 ;
    std::vector<long> tmp_v_index;
    for(int i = range_low_x ; i <= range_high_x; i++ )
    {
        for(int j=range_low_mom ; j < range_high_mom ; j++)
        {
            for(long k=0 ; k < conts.at(i).at(j).size() ; k++)
            {
                long in = conts.at(i).at(j).at(k);
                if( abs(db_px[in]-q_px) > float_value || abs(db_py[in]-q_py) > float_value || abs(db_pz[in]-q_pz) > float_value ) continue;
                tmp_v_index.push_back(in);
            }
        }
    }
    if(tmp_v_index.size()==0) return -1;
    else {
        int ran = rand()%tmp_v_index.size();
        long index = tmp_v_index.at(ran);
        //en_scale = tmp_mom/sqrt(db_px[index]*db_px[index] + db_py[index]*db_py[index] + db_pz[index]*db_pz[index]) ;
        en_scale = 1;//FIXME check 
        b_x=db_x[index]; b_px=db_px[index]; b_py=db_py[index]; b_pz=db_pz[index];
        //std::cout<<"tmp_v_index_size="<<tmp_v_index.size()<<",ran-"<<ran<<",index="<<index<<",en_scale="<<en_scale<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<std::endl;
        return index; 
    }
}

long FrozenShowerPlugin::Search_v1(const float& low_x, const float& low_theta, const float& low_phi, const float& low_mom, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector< std::vector< std::vector<long> > > > > & conts, float& en_scale, float& b_x, float& b_px, float& b_py, float& b_pz, TH1D* h_res)
{ 
    float tmp_mom = sqrt(q_px*q_px + q_py*q_py  + q_pz*q_pz);
    float tmp_theta = acos(q_pz/tmp_mom)*180/PI ;
    float tmp_phi = acos(q_px/sqrt(q_px*q_px + q_py*q_py) )*180/PI ;
    if(q_py<0) tmp_phi = -tmp_phi;
    int flow_x = 0 ;
    //int flow_theta = 1 ;
    //int flow_phi   = 1 ;
    int flow_theta = 0 ;
    int flow_phi   = 0 ;
    int flow_mom = int(0.1*tmp_mom) ;
    //int flow_mom = int(0.3*tmp_mom) ;
    //int flow_mom = int(0.2*tmp_mom) ;
    //int flow_mom = int(0.9*tmp_mom) ;
    
    //float float_value = 0.1*tmp_mom;
    //int flow_mom = int(0.1*tmp_mom) ;
    //float float_value = 0.02*tmp_mom;
    //int flow_x = 0 ;
    //int flow_mom = int(0.03*tmp_mom) ;
    int range_low_x  = int(q_x-low_x-flow_x) > 0            ? int(q_x-low_x-flow_x) :            0   ;
    int range_high_x = int(q_x-low_x+flow_x) < conts.size() ? int(q_x-low_x+flow_x) : conts.size()-1 ;
    int range_low_theta  = int(tmp_theta-low_theta-flow_theta) > 0                  ? int(tmp_theta-low_theta-flow_theta) :            0         ;
    int range_high_theta = int(tmp_theta-low_theta+flow_theta) < conts.at(0).size() ? int(tmp_theta-low_theta+flow_theta) : conts.at(0).size()-1 ;
    int range_low_phi  = int(tmp_phi-low_phi-flow_phi) > 0                  ? int(tmp_phi-low_phi-flow_phi) :            0                     ;
    int range_high_phi = int(tmp_phi-low_phi+flow_phi) < conts.at(0).at(0).size() ? int(tmp_phi-low_phi+flow_phi) : conts.at(0).at(0).size()-1 ;
    int range_low_mom  = int(tmp_mom-low_mom-flow_mom) > 0                              ? int(tmp_mom-low_mom-flow_mom) :                  0               ;
    int range_high_mom = int(tmp_mom-low_mom+flow_mom) < conts.at(0).at(0).at(0).size() ? int(tmp_mom-low_mom+flow_mom) : conts.at(0).at(0).at(0).size()-1 ;
    std::vector<long> tmp_v_index;
    //std::cout<<"x="<<q_x<<",px="<<q_px<<",py="<<q_py<<",pz="<<q_pz<<",theta="<<tmp_theta<<",phi="<<tmp_phi<<",low_x="<<range_low_x<<",high_x="<<range_high_x<<",low_theta="<<range_low_theta<<",high_theta="<<range_high_theta<<",low_phi="<<range_low_phi<<",high_phi="<<range_high_phi<<",low_mom="<<range_low_mom<<",high_mom="<<range_high_mom<<std::endl;
    for(int i = range_low_x ; i <= range_high_x; i++ )
    {
        for(int ii = range_low_theta ; ii <= range_high_theta; ii++ )
        {
            for(int iii = range_low_phi ; iii <= range_high_phi; iii++ )
            {
                for(int j=range_low_mom ; j < range_high_mom ; j++)
                {
                    for(long k=0 ; k < conts.at(i).at(ii).at(iii).at(j).size() ; k++)
                    {
                        long in = conts.at(i).at(ii).at(iii).at(j).at(k);
                        tmp_v_index.push_back(in);
                    }
                }
            }
        }
    }
    //std::cout<<"tmp_v_index_size="<<tmp_v_index.size()<<std::endl;
    if(tmp_v_index.size()==0) return -1;
    else {
        int ran = rand()%tmp_v_index.size();
        long index = tmp_v_index.at(ran);
        en_scale = tmp_mom/sqrt(db_px[index]*db_px[index] + db_py[index]*db_py[index] + db_pz[index]*db_pz[index]) ;
        //float rang = gaussrand()/10.0;// res 0.1
        //float rang = tmp_mom > 500 ? gaussrand()/5.0 : gaussrand()/10.0;// res 0.1
        int tmp_bin = h_res->FindBin(double(tmp_mom));
        float tmp_res = h_res->GetBinContent(tmp_bin);
        //float rang = tmp_mom > 500 ? CLHEP::RandGauss::shoot(0, 0.2) : CLHEP::RandGauss::shoot(0, 0.1);// res 0.1
        //float rang = CLHEP::RandGauss::shoot(0, tmp_res);// res 0.1
        //float rang = tmp_mom < 800 ? CLHEP::RandGauss::shoot(0, 0.1) : CLHEP::RandGauss::shoot(0, 0.3);// res 0.1
        //float rang = tmp_mom < 800 ? CLHEP::RandGauss::shoot(0, 0.15) : CLHEP::RandGauss::shoot(0, 0.05);// res 0.1
        //float rang = CLHEP::RandGauss::shoot(0, 0.1);//nominal
        float rang = 0;// check 
        //std::cout<<"rang="<<rang<<std::endl;
        en_scale = en_scale*(1+rang); 
        //en_scale = 1;//FIXME check 
        b_x=db_x[index]; b_px=db_px[index]; b_py=db_py[index]; b_pz=db_pz[index];
        //std::cout<<"tmp_v_index_size="<<tmp_v_index.size()<<",ran-"<<ran<<",index="<<index<<",en_scale="<<en_scale<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<std::endl;
        return index; 
    }
}

void FrozenShowerPlugin::processTxt(const std::string& mc_info, float* & array_0, float*  & array_1, float* & array_2, float* & array_3, bool isBarrel, std::vector< std::vector< std::vector<long> > > & v_v_vec) {
  int low = 0;
  int high = 0;
  if(isBarrel){
      low  = m_x_min;
      high = m_x_max;
  }
  else{
      low  = abs_endcap_z_min;
      high = abs_endcap_z_max;
  }

  for(int i=low; i<=high; i++)//1 mm
  {
      std::vector< std::vector<long> > tmp_x;
      for(int j=m_Mom_min; j<=m_Mom_max; j++)// 1 MeV
      {
          std::vector<long> tmp_Mom;
          tmp_x.push_back(tmp_Mom);
      }
      v_v_vec.push_back(tmp_x);
  }


  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  long N_total = 0;
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  while(getline(m_f,sline)){
      N_total ++ ;
  }
  m_f.close();
  std::cout << "total entry ="<<N_total<<" in "<<mc_info<< std::endl;
  array_0 = new float[N_total];
  array_1 = new float[N_total];
  array_2 = new float[N_total];
  array_3 = new float[N_total];
  int tmp_N = 0;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      array_0[tmp_N] = ::atof(s1.c_str());
      array_1[tmp_N] = ::atof(s2.c_str());
      array_2[tmp_N] = ::atof(s3.c_str());
      array_3[tmp_N] = ::atof(s4.c_str());
      float tmp_mom =sqrt( array_1[tmp_N]*array_1[tmp_N] + array_2[tmp_N]*array_2[tmp_N] + array_3[tmp_N]*array_3[tmp_N] );
      if( isBarrel && array_0 [tmp_N] > m_x_min  && array_0 [tmp_N] < m_x_max && tmp_mom > m_Mom_min && tmp_mom < m_Mom_max  ) v_v_vec.at(int(::atof(s1.c_str())-m_x_min)).at(int(tmp_mom-m_Mom_min)).push_back(tmp_N);
      else if( isBarrel==false && array_0 [tmp_N] > abs_endcap_z_min  && array_0 [tmp_N] < abs_endcap_z_max && tmp_mom > m_Mom_min && tmp_mom < m_Mom_max  ) v_v_vec.at(int(::atof(s1.c_str())-abs_endcap_z_min)).at(int(tmp_mom-m_Mom_min)).push_back(tmp_N);
      tmp_N ++ ;
  }
}

void FrozenShowerPlugin::createMomHist(TH1D* h1, float* px, float* py, float* pz, long n_tot) {
    std::cout<<"n_tot size="<<n_tot<<std::endl;
    for(unsigned int i=0; i< n_tot; i++){
        float tmp_mom =sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] );
        h1->Fill(tmp_mom);
    }
}
void FrozenShowerPlugin::printHist(const TH1D* h1 ) {
    for(unsigned int i=1; i <= h1->GetNbinsX(); i++){
        std::cout<<"bin="<<i<<",val="<<h1->GetBinContent(i)<<std::endl;
    }
}

void FrozenShowerPlugin::createResHist(const TH1D* h1, TH1D* h_res, float res ) {
    long tmp_min = 1000000 ;
    int tmp_index = -1;
    for(unsigned int i=1; i <= h1->GetNbinsX(); i++){
        if(h1->GetBinContent(i) > 100){
            if(h1->GetBinContent(i)<tmp_min) { tmp_min = h1->GetBinContent(i); tmp_index=i;}
        }
    }
    if(tmp_index != -1) std::cout<<"min="<<tmp_min<<",bin="<<tmp_index<<std::endl;
    else throw("Error, too small sps!");
    float tmp_c = sqrt(tmp_min)*res;
    for (unsigned int i=1; i <= h1->GetNbinsX(); i++){
        h_res->SetBinContent(i, tmp_c/sqrt(h1->GetBinContent(i)));
    }   
}

void FrozenShowerPlugin::processTxt(const std::string& mc_info, float* & array_0, float*  & array_1, float* & array_2, float* & array_3, bool isBarrel, std::vector< std::vector< std::vector< std::vector< std::vector<long> > > > > & v_v_vec, long & n_tot) {
  int low = 0;
  int high = 0;
  if(isBarrel){
      low  = m_x_min;
      high = m_x_max;
  }
  else{
      low  = abs_endcap_z_min;
      high = abs_endcap_z_max;
  }

  for(int i=low; i<=high; i++)//1 mm
  {
      std::vector< std::vector< std::vector< std::vector<long> > > > tmp_x;
      for(int j=m_theta_min; j<=m_theta_max; j++)// 1 MeV
      {
          std::vector< std::vector< std::vector<long> > > tmp_theta;
          for(int k=m_phi_min; k<=m_phi_max; k++)// 1 MeV
          {
              std::vector< std::vector<long> > tmp_phi;
              for(int z=m_Mom_min; z<=m_Mom_max; z++)// 1 MeV
              {
                  std::vector<long> tmp_Mom;
                  tmp_phi.push_back(tmp_Mom);
              }
              tmp_theta.push_back(tmp_phi);
          }
          tmp_x.push_back(tmp_theta);
      }
      v_v_vec.push_back(tmp_x);
  }


  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  long N_total = 0;
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  while(getline(m_f,sline)){
      N_total ++ ;
  }
  m_f.close();
  std::cout << "total entry ="<<N_total<<" in "<<mc_info<< std::endl;
  n_tot = N_total;
  array_0 = new float[N_total];
  array_1 = new float[N_total];
  array_2 = new float[N_total];
  array_3 = new float[N_total];
  int tmp_N = 0;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      array_0[tmp_N] = ::atof(s1.c_str());
      array_1[tmp_N] = ::atof(s2.c_str());
      array_2[tmp_N] = ::atof(s3.c_str());
      array_3[tmp_N] = ::atof(s4.c_str());
      float tmp_mom =sqrt( array_1[tmp_N]*array_1[tmp_N] + array_2[tmp_N]*array_2[tmp_N] + array_3[tmp_N]*array_3[tmp_N] );
      float tmp_theta = acos( array_3[tmp_N]/tmp_mom ) * 180 /PI ;
      float tmp_phi = acos( array_1[tmp_N]/sqrt( array_1[tmp_N]*array_1[tmp_N] + array_2[tmp_N]*array_2[tmp_N] ) ) * 180 /PI ;
      if(array_2[tmp_N] < 0) tmp_phi = -tmp_phi;
      if( isBarrel && array_0 [tmp_N] > m_x_min  && array_0 [tmp_N] < m_x_max && tmp_theta > m_theta_min && tmp_theta < m_theta_max && tmp_phi > m_phi_min && tmp_phi < m_phi_max && tmp_mom > m_Mom_min && tmp_mom < m_Mom_max  ) v_v_vec.at(int(::atof(s1.c_str())-m_x_min)).at(int(tmp_theta-m_theta_min)).at(int(tmp_phi-m_phi_min)).at(int(tmp_mom-m_Mom_min)).push_back(tmp_N);
      //else if( isBarrel==false && array_0 [tmp_N] > abs_endcap_z_min  && array_0 [tmp_N] < abs_endcap_z_max && tmp_mom > m_Mom_min && tmp_mom < m_Mom_max  ) v_v_vec.at(int(::atof(s1.c_str())-abs_endcap_z_min)).at(int(tmp_mom-m_Mom_min)).push_back(tmp_N);
      tmp_N ++ ;
  }
}

void FrozenShowerPlugin::processMapping( map<int, TChain*>& chain_map, const std::string& x_s, const std::string& file) {

     TChain* tmp_tree = new TChain("evt");
     tmp_tree->Add(file.c_str());
     tmp_tree->SetBranchAddress("m_ECAL_Hit_x", &m_ECAL_Hit_x       );
     tmp_tree->SetBranchAddress("m_ECAL_Hit_y", &m_ECAL_Hit_y       );
     tmp_tree->SetBranchAddress("m_ECAL_Hit_z", &m_ECAL_Hit_z       );
     tmp_tree->SetBranchAddress("m_ECAL_Hit_E", &m_ECAL_Hit_E       );
     tmp_tree->SetBranchAddress("m_Tag"       , &m_Tag              );
     tmp_tree->SetBranchAddress("m_ECAL_ID_I" , &m_ECAL_ID_I        );
     tmp_tree->SetBranchAddress("m_ECAL_ID_K" , &m_ECAL_ID_K        );
     if(x_s=="e-")
     {
         chain_map.insert(pair<int,  TChain*>(11 ,tmp_tree));    
         std::cout<<"total_entry_em="<<tmp_tree->GetEntries()<<std::endl;
     }
     else if(x_s=="e+") 
     {
         chain_map.insert(pair<int,  TChain*>(-11 ,tmp_tree));    
         std::cout<<"total_entry_ep="<<tmp_tree->GetEntries()<<std::endl;
     }
     else{
         std::cout<<"wrong name, break"<<std::endl;
         throw ;
     }
     std::cout << "Done Set library..."<<std::endl;
}

