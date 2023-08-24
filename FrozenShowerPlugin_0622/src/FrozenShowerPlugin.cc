// ********************************************************
// *            Mokka    http://mokka.in2p3.fr            *
// *    -- A Detailed Geant4 Simulation for the ILC --    *
// ********************************************************
//
// $Id: FrozenShowerPlugin.cc,v 1.3 2008/02/08 14:01:24 adrian Exp $
// $Name: mokka-06-06 $
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




// hit position from cell center
// using faiss
// library for electron and positron, pid_bin =0 for 11, pid_bin=1 for -11
// Using library from start points
// E 900 bins for 100MeV to 1000MeV, phi 50 bins for -25 to 25, theta 40 bins for |theta-90|<40.
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
  doSaveSP    = false;
  checkTime   = true;
  saveNotFind = false;

  theEncoder= new Encoder32();
  width_M = 2350*2/5;
  width_J = width_M/90;

  x_min = 1850;
  x_max = 2000;//mm
  theta_min = 0;
  theta_max = 180;
  phi_min = 0;
  phi_max = 360;
  Mom_min = 100;//MeV
  Mom_max = 1000;//MeV

  y_boundry = 500;//mm 765
  z_boundry = 2200;//mm 2350
  
  outCut = 0.02;//2%
  //std::istringstream string_rad((*Control::globalModelParameters)["tracker_region_rmax"]);
  //std::istringstream string_z((*Control::globalModelParameters)["tracker_region_zmax"]);

  //string_rad >> tracking_radius_max;
  //string_z >> tracking_z_max;
  if(doSaveSP) std::cout << "Using save start point model" << std::endl;
  else{

  for(int i=x_min; i<=x_max; i++)
  {
      std::vector<long> tmp;
      m_em_x_vec.push_back(tmp);
      m_ep_x_vec.push_back(tmp);
  }


  std::ifstream m_f_em;
  std::ifstream m_f_ep;
  std::string mc_info_em = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_addNoFind//lib_e-.txt";
  std::string mc_info_ep = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_addNoFind//lib_e+.txt";
  //std::string mc_info_em = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0522//lib_e-.txt";
  //std::string mc_info_ep = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0522//lib_e+.txt";
  //std::string mc_info_em = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0529_smear//lib_e-.txt";
  //std::string mc_info_ep = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0529_smear//lib_e+.txt";
  m_f_em.open(mc_info_em.c_str());
  if(m_f_em.fail()) { std::cout << "error, can't open "<<mc_info_em<< std::endl; throw ;}
  m_em_total = 0;
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  while(getline(m_f_em,sline)){
      m_em_total ++ ;
  }
  m_f_em.close();
  std::cout << "total entry ="<<m_em_total<<" in "<<mc_info_em<< std::endl;
  m_em_x  = new float[m_em_total];
  m_em_px = new float[m_em_total];
  m_em_py = new float[m_em_total];
  m_em_pz = new float[m_em_total];
  int tmp_N = 0;
  m_f_em.open(mc_info_em.c_str());
  if(m_f_em.fail()) { std::cout << "error, can't open "<<mc_info_em<< std::endl; throw ;}
  while(getline(m_f_em,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      m_em_x [tmp_N] = ::atof(s1.c_str());
      m_em_px[tmp_N] = ::atof(s2.c_str());
      m_em_py[tmp_N] = ::atof(s3.c_str());
      m_em_pz[tmp_N] = ::atof(s4.c_str());
      if     (m_em_x [tmp_N]<1900 && m_em_py[tmp_N]>0)   m_em_cont_x1900L_pyP.push_back(tmp_N);
      else if(m_em_x [tmp_N]<1900 && m_em_py[tmp_N]<=0)  m_em_cont_x1900L_pyM.push_back(tmp_N);
      else if(m_em_x [tmp_N]>=1900 && m_em_py[tmp_N]>0)  m_em_cont_x1900H_pyP.push_back(tmp_N);
      else if(m_em_x [tmp_N]>=1900 && m_em_py[tmp_N]<=0) m_em_cont_x1900H_pyM.push_back(tmp_N);
      if(::atof(s1.c_str()) > x_min  && ::atof(s1.c_str()) < x_max )  m_em_x_vec.at(int(::atof(s1.c_str())-x_min)).push_back(tmp_N);
      tmp_N ++ ;
  }
  //////////////////
  m_f_ep.open(mc_info_ep.c_str());
  if(m_f_ep.fail()) { std::cout << "error, can't open "<<mc_info_ep<< std::endl; throw ;}
  m_ep_total = 0;
  while(getline(m_f_ep,sline)){
      m_ep_total ++ ;
  }
  m_f_ep.close();
  std::cout << "total entry ="<<m_ep_total<<" in "<<mc_info_ep<< std::endl;
  m_ep_x  = new float[m_ep_total];
  m_ep_px = new float[m_ep_total];
  m_ep_py = new float[m_ep_total];
  m_ep_pz = new float[m_ep_total];
  tmp_N = 0;
  m_f_ep.open(mc_info_ep.c_str());
  if(m_f_ep.fail()) { std::cout << "error, can't open "<<mc_info_ep<< std::endl; throw ;}
  while(getline(m_f_ep,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      m_ep_x [tmp_N] = ::atof(s1.c_str());
      m_ep_px[tmp_N] = ::atof(s2.c_str());
      m_ep_py[tmp_N] = ::atof(s3.c_str());
      m_ep_pz[tmp_N] = ::atof(s4.c_str());
      if     (m_ep_x [tmp_N]<1900  && m_ep_py[tmp_N]>0)  m_ep_cont_x1900L_pyP.push_back(tmp_N);
      else if(m_ep_x [tmp_N]<1900  && m_ep_py[tmp_N]<=0) m_ep_cont_x1900L_pyM.push_back(tmp_N);
      else if(m_ep_x [tmp_N]>=1900 && m_ep_py[tmp_N]>0)  m_ep_cont_x1900H_pyP.push_back(tmp_N);
      else if(m_ep_x [tmp_N]>=1900 && m_ep_py[tmp_N]<=0) m_ep_cont_x1900H_pyM.push_back(tmp_N);
      if(::atof(s1.c_str()) > x_min  && ::atof(s1.c_str()) < x_max )  m_ep_x_vec.at(int(::atof(s1.c_str())-x_min)).push_back(tmp_N);
      tmp_N ++ ;
  }

  /////// 
  std::cout << "Set library..."<<std::endl;
  std::vector<std::string> v_string;
  v_string.push_back("e-");
  v_string.push_back("e+");
  for(unsigned int i =0 ; i < v_string.size(); i++)
  {
     TChain* tmp_tree = new TChain("evt");
     std::string x_s = v_string.at(i);
     std::string file = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_addNoFind/lib_"+x_s+".root";
     //std::string file = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0522/lib_"+x_s+".root";
     //std::string file = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0529_smear/lib_"+x_s+".root";
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
         x_chain_map.insert(pair<int,  TChain*>(11 ,tmp_tree));    
         total_entry_em = tmp_tree->GetEntries();
         std::cout<<"total_entry_em="<<total_entry_em<<std::endl;
     }
     else if(x_s=="e+") 
     {
         x_chain_map.insert(pair<int,  TChain*>(-11 ,tmp_tree));    
         total_entry_ep = tmp_tree->GetEntries();
         std::cout<<"total_entry_ep="<<total_entry_ep<<std::endl;
     }
     else{
         std::cout<<"wrong name, break"<<std::endl;
         throw ;
     }
  }
  std::cout << "Done Set library..."<<std::endl;

  }//else
  if(checkTime)
  {
      std::cout << "Check time..."<<std::endl;
      std::vector<float>* hits_x;
      std::vector<float>* hits_y;
      std::vector<float>* hits_z;
      std::vector<float>* hits_E;
      int  tag;
      std::vector<int>* hits_I;
      std::vector<int>* hits_K;
      for(unsigned int i = 0; i< 100000; i++)
      {
          x_chain_map[11]->GetEntry(i);
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
          x_chain_map[-11]->GetEntry(i);
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

  int Depth[29] = {1850, 1857, 1860, 1868, 1871, 1878, 1881, 1889, 1892, 1899, 1902, 1910, 1913, 1920, 1923, 1931, 1934, 1941, 1944, 1952, 1957, 1967, 1972, 1981, 1986, 1996, 2001, 2011, 2016};
  for(unsigned int i=0; i<29; i++) m_layer.push_back(Depth[i]);

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

  //x_chain_map[1870]->GetEntry(int(10*90*400+10*400+200));
  //std::cout<<"m_ECAL_ID_I size="<<m_ECAL_ID_I->size()<<std::endl;
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
int pid = 0;
double x = 0;
double y = 0;
double z = 0;
double px = 0;
double py = 0;
double pz = 0;
int a = makeHitCollection(step, pid, x, y, z, px, py, pz);
m_return.push_back(a);
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
  clock_t t0, t1;
  //t0 = clock();
  G4Track* theTrack  =  step->GetTrack(); 
  const G4ParticleDefinition* def = theTrack->GetDefinition();
  G4int pdgCode=def->GetPDGEncoding();
  _pid = pdgCode ;
  
  const G4ThreeVector &postPos = step->GetPostStepPoint()->GetPosition();
  const G4ThreeVector &postMom = step->GetPostStepPoint()->GetMomentum();
  double radius_sqrd = postPos[0]*postPos[0] + postPos[1]*postPos[1];
  float Mom  = sqrt(postMom[0]*postMom[0] + postMom[1]*postMom[1] + postMom[2]*postMom[2]);
  _x = postPos[0];
  _y = postPos[1];
  _z = postPos[2];
  _px = postMom[0];
  _py = postMom[1];
  _pz = postMom[2];
  if( abs(pdgCode) != 11 ) return 0;
  if( radius_sqrd < x_min*x_min || abs(postPos[2]) > z_boundry || Mom < Mom_min || Mom > Mom_max ) return 2;

  G4double time = step->GetTrack()->GetGlobalTime() ;
  const float PI = acos(-1);

  
  float phi = getPhi(postPos[0], postPos[1]);
  float phi_p  = getPhi(postMom[0], postMom[1]);
  int part = partition(phi);// phi should < 360 and > 0
  //int Stave = part>=3 ? (part-3) : (part+5); 
  int Stave = part>=3 ? (part-2) : (part+6); 
  //std::cout <<"part="<<part<<"phi="<<phi<<",x="<<postPos[0]<<",y="<<postPos[1]<<  std::endl;
  float rotated = (part-1)*45;
  float tmp_phi = abs(phi-rotated) < 22.6 ? abs(phi-rotated) : abs(360-(phi-rotated));
  //if(tmp_phi>15) return 2;//gap
  float r = sqrt(postPos[0]*postPos[0] + postPos[1]*postPos[1]);
  float pt = sqrt(postMom[0]*postMom[0] + postMom[1]*postMom[1]);
  float new_x = r*cos((phi-rotated)*PI/180);
  float new_y = r*sin((phi-rotated)*PI/180);
  float new_px = pt*cos((phi_p-rotated)*PI/180);
  float new_py = pt*sin((phi_p-rotated)*PI/180);
  
  float cos_theta = postMom[2] / Mom;
  float theta = acos(cos_theta)*180/PI;//0-180
  //if(new_px < 0 || new_x<1870 || new_x>1997 ) return 3; //select  
  if(new_px < 0 || new_x<1851 || new_x>2000 ) return 3; //select  
  if( abs(new_y) > y_boundry ) return 3; //from not found start points, don't consider edge now  
  if(doSaveSP)
  {
      //m_point_x.push_back(new_x);
      //m_point_y.push_back(new_y);
      m_point_x.push_back(postPos[0]);
      m_point_y.push_back(postPos[1]);
      m_point_z.push_back(postPos[2]);
      //m_mom_x  .push_back(new_px);
      //m_mom_y  .push_back(new_py);
      m_mom_x  .push_back(postMom[0]);
      m_mom_y  .push_back(postMom[1]);
      m_mom_z  .push_back(postMom[2]);
      m_pid    .push_back(pdgCode   );
      theTrack->SetTrackStatus(fStopAndKill);
      return 111;
  }
  /*
  int bin_pid = (pdgCode == 11) ? 0 : 1 ;  
  float cos_phi = new_px / pt;
  int cos_phi_1 = int(acos(cos_phi)*180/PI);
  if(cos_phi_1>=25) return 4;
  int bin_phi = (new_py > 0) ? cos_phi_1 : 25 + cos_phi_1;//
  int bin_theta = int(abs(theta-90));
  if(bin_theta>=40) return 5;
  int bin_E = int(Mom-100);//already be MeV
  if(bin_E>=900) return 6;
  int index = bin_pid*50*40*900 + bin_phi*40*900 + bin_theta*900 + bin_E;
  */
  long nns = -1; 
  //t0 = clock();
  if(pdgCode==11)
  {
      //if     (new_x < 1900 && new_py > 0)    nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_em_x, m_em_px, m_em_py , m_em_pz, m_em_cont_x1900L_pyP ) ;
      //else if(new_x < 1900  && new_py <= 0)  nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_em_x, m_em_px, m_em_py , m_em_pz, m_em_cont_x1900L_pyM ) ;
      //else if(new_x >= 1900 && new_py > 0 )  nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_em_x, m_em_px, m_em_py , m_em_pz, m_em_cont_x1900H_pyP ) ;
      //else if(new_x >= 1900 && new_py <= 0 ) nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_em_x, m_em_px, m_em_py , m_em_pz, m_em_cont_x1900H_pyM ) ;
      nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_em_x, m_em_px, m_em_py , m_em_pz, m_em_x_vec) ;
      //nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_em_x, m_em_px, m_em_py , m_em_pz,m_em_total) ;
      if(nns != -1 )
      {
          if(checkTime==false){x_chain_map[11]->GetEntry(nns);
                               //std::cout<<"em nns="<<nns<<std::endl;
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
              m_noFind_point_x.push_back(new_x);
              m_noFind_point_y.push_back(new_y);
              m_noFind_point_z.push_back(postPos[2]);
              m_noFind_mom_x  .push_back(new_px);
              m_noFind_mom_y  .push_back(new_py);
              m_noFind_mom_z  .push_back(postMom[2]);
              m_noFind_pid    .push_back(pdgCode   );
          }
      return 4;
      }
  }
  else if(pdgCode==-11)
  {
      //if     (new_x < 1900 && new_py > 0)    nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_ep_x, m_ep_px, m_ep_py , m_ep_pz, m_ep_cont_x1900L_pyP ) ;
      //else if(new_x < 1900  && new_py <= 0)  nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_ep_x, m_ep_px, m_ep_py , m_ep_pz, m_ep_cont_x1900L_pyM ) ;
      //else if(new_x >= 1900 && new_py > 0 )  nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_ep_x, m_ep_px, m_ep_py , m_ep_pz, m_ep_cont_x1900H_pyP ) ;
      //else if(new_x >= 1900 && new_py <= 0 ) nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_ep_x, m_ep_px, m_ep_py , m_ep_pz, m_ep_cont_x1900H_pyM ) ;
      nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_ep_x, m_ep_px, m_ep_py , m_ep_pz, m_ep_x_vec ) ;
      //nns = Search(new_x, new_px, new_py, abs(postMom[2]), m_ep_x, m_ep_px, m_ep_py , m_ep_pz, m_ep_total ) ;
      if(nns != -1 ) 
      {
          if(checkTime==false){x_chain_map[-11]->GetEntry(nns);
                               //std::cout<<"ep nns="<<nns<<std::endl;
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
              m_noFind_point_x.push_back(new_x);
              m_noFind_point_y.push_back(new_y);
              m_noFind_point_z.push_back(postPos[2]);
              m_noFind_mom_x  .push_back(new_px);
              m_noFind_mom_y  .push_back(new_py);
              m_noFind_mom_z  .push_back(postMom[2]);
              m_noFind_pid    .push_back(pdgCode   );
          }
          return 5;
      }
  }
  //t1 = clock();
  //std::cout << "Hi search, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  
  m_pass0 ++ ;
  if(m_Tag == 0) return 6;
  m_pass1 ++ ;
  //t1 = clock();
  //std::cout << "Hi 1, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  //t0 = clock();
  //# time consuming part !! ####
  //map<int, TChain*>::iterator it = x_chain_map.find(pdgCode);
  //if( it != x_chain_map.end()){it->second->GetEntry(index);}
  //else {return 8;}
  
  //std::cout<<"new_x="<<new_x<<",int(new_x)="<<int(new_x)<<",pid="<<_pid<<",x="<<_x<<",y="<<_y<<",z"<<_z<<",_px"<<_px<<",py="<<_py<<",pz="<<_pz<<",bin_theta="<<bin_theta<<",bin_phi="<<bin_phi<<",index="<<index<<",m_Tag->size()="<<m_Tag->size()<<std::endl; 
  //##############
  /*
  t1 = clock();
  std::cout << "Hi map, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  t0 = clock();
  */
  float tot_E = 0;
  float out_boundry_E = 0;
  //++++++++++++++++++++++++++++++++++++
  //std::cout <<"After index, m_ECAL_Hit_y size="<<m_ECAL_Hit_y->size()<<",m_Tag->size()="<<m_Tag->size()<<",lib_index="<<lib_index<<std::endl;
  //########################################
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SEcalSD02* theSD = dynamic_cast<SEcalSD02*> ( SDman->FindSensitiveDetector("EcalBarrelSilicon") ); 
  if(!theSD) {std::cout<<"NotFoundSD EcalBarrelSilicon"<<std::endl; return 11;}
  HitsCollection *CalCollection;
  CalCollection = theSD->NormalCalCollection;
  G4int n_hit = CalCollection->entries();
  G4int PDG = pdgCode;
  G4Step* astep = const_cast<G4Step*>(step);
  G4int PID =  Control::GetControl()->GetPIDForCalHit(astep);//for found primary MC particle

  //#####################################
  /*
  t1 = clock();
  std::cout << "Hi 2, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  t0 = clock();
  */
  G4int theSDPiece = 0;
  theSDPiece = ECALBARREL;
  for(unsigned int i=0; i< m_ECAL_Hit_y->size(); i++)
  {
      float Hit_x     = m_ECAL_Hit_x->at(i) + new_x ;
      int l_index = -1;
      //l_index = m_ECAL_ID_K->at(i);
      l_index = m_ECAL_ID_K->at(i) + 2;//due to get from K-1
      if(l_index==-1) continue;
      //std::cout<<"l_index="<<l_index<<",Hit_x="<<Hit_x<<std::endl;
      float Hit_y     = m_ECAL_Hit_y->at(i) ;
      float Hit_z     = m_ECAL_Hit_z->at(i) ;
      // parity with x-z plane 
      //float new_Hit_y = new_py > 0 ? Hit_y : -Hit_y ;
      float new_Hit_y = Hit_y ;
      // parity with x-y plane 
      float new_Hit_z = postMom[2] > 0 ? Hit_z : -Hit_z ;
      //moving
      new_Hit_y = new_Hit_y + new_y;//from center to new y
      int I_index = m_ECAL_ID_I->at(i) - int(new_y/10);// new I
      new_Hit_z = new_Hit_z + postPos[2];
      tot_E += m_ECAL_Hit_E->at(i);
      //if(abs(new_Hit_y) > y_boundry || abs(new_Hit_z) > z_boundry) out_boundry_E += m_ECAL_Hit_E->at(lib_index).at(i);
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
                  theCellCenter = theSD->GetCellCenter(theSDPiece, Stave, tmp_M, tmp_I, tmp_J, l_index); 
                  dis = theCellCenter - thePosition ;
                  if ( (dis.mag() < min_dis) && (dis.mag() < theSD->CellDim.mag()*sqrt(2.))) 
                  {
                      secd_min_dis = min_dis;
                      min_dis = dis.mag();
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
      cell_ids theIndex =  theEncoder->encode(Stave, final_M, final_I, final_J, l_index, 0);
      G4ThreeVector finalCellCenter = theSD->GetCellCenter(theSDPiece, Stave, final_M, final_I, final_J, l_index);
      if(found_id == false)
      {
          finalCellCenter =  thePosition ;
      }
      //std::cout<<"min_dis="<<min_dis<<",secd mis="<<secd_min_dis<<",found_id="<<found_id<<",CellDim.mag()="<<theSD->CellDim.mag()<<",x="<<finalCellCenter(0)<<",y="<<finalCellCenter(1)<<",z="<<finalCellCenter(2)<<",real x="<<real_Hit_x<<", real y="<<real_Hit_y<<", real z="<<real_Hit_z<< std::endl;
      float edep = m_ECAL_Hit_E->at(i)*1000 ;//to MeV
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
				  Stave      ,
				  final_M,
				  final_I,
				  final_J,
				  l_index,
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
  theTrack->SetTrackStatus(fStopAndKill);
  //t1 = clock();
  //std::cout << "Hi end, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  return 12;      
}

float FrozenShowerPlugin::getPhi(float x, float y)
{
    if     (x==0 && y>0) return 90;
    else if(x==0 && y<0) return 270;
    else if(x==0 && y==0) return 0;
    float phi = atan(y/x)*180/acos(-1);
    if                 (x<0) phi = phi + 180;
    else if     (x>0 && y<0) phi = phi + 360;
    return phi;
}

int FrozenShowerPlugin::partition(float phi) // phi should < 360 and > 0
{   if(phi<0) throw "Wrong input phi!";
    if(phi<=22.5 || phi > (360-22.5)) return 1;
    else if(22.5       <phi && phi<=(22.5+1*45)) return 2;
    else if((22.5+1*45)<phi && phi<=(22.5+2*45)) return 3;
    else if((22.5+2*45)<phi && phi<=(22.5+3*45)) return 4;
    else if((22.5+3*45)<phi && phi<=(22.5+4*45)) return 5;
    else if((22.5+4*45)<phi && phi<=(22.5+5*45)) return 6;
    else if((22.5+5*45)<phi && phi<=(22.5+6*45)) return 7;
    else if((22.5+6*45)<phi && phi<=(22.5+7*45)) return 8;
    else{std::cout<<"something wrong"<<std::endl;return -1;}
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

/*
long FrozenShowerPlugin::Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector<long> > & conts) 
{  
    long index=-1;
    float min_dist = 5;
    int flow = 1 ;
    int range_low  = int(q_x-x_min-flow) > 0            ? int(q_x-x_min-flow) :            0   ;
    int range_high = int(q_x-x_min+flow) < conts.size() ? int(q_x-x_min+flow) : conts.size()-1 ;
    for(int i = range_low ; i <= range_high; i++ )
    {
        //std::cout<<"in_size="<<conts.at(i).size()<<std::endl;
        for(long k=0 ; k < conts.at(i).size() ; k++)
        {
            long in = conts.at(i).at(k);
            if(abs(db_x[in]-q_x) > min_dist || abs(db_px[in]-q_px) > min_dist || abs(db_py[in]-q_py) > min_dist || abs(db_pz[in]-q_pz) > min_dist ) continue;
            float dist =  abs(db_x[in]-q_x) + abs(db_px[in]-q_px) + abs(db_py[in]-q_py) + abs(db_pz[in]-q_pz) ;
            if( dist < min_dist ) {min_dist = dist; index = i;}
            if( min_dist < 1) break;
        }
    }
    //std::cout<<"conts_size="<<conts.size()<<",index="<<index<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<",dist="<< min_dist <<std::endl;
    return index; 
}
*/


long FrozenShowerPlugin::Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector<long> > & conts) 
{ 
    long index=-1;
    float min_dist = 10;
    int flow = 1 ;
    int range_low  = int(q_x-x_min-flow) > 0            ? int(q_x-x_min-flow) :            0   ;
    int range_high = int(q_x-x_min+flow) < conts.size() ? int(q_x-x_min+flow) : conts.size()-1 ;
    for(int i = range_low ; i <= range_high; i++ )
    {
        //std::cout<<"in_size="<<conts.at(i).size()<<std::endl;
        for(long k=0 ; k < conts.at(i).size() ; k++)
        {
            long in = conts.at(i).at(k);
            float pre_cut = sqrt(min_dist);
            if(abs(db_x[in]-q_x) > pre_cut || abs(db_px[in]-q_px) > pre_cut || abs(db_py[in]-q_py) > pre_cut || abs(db_pz[in]-q_pz) > pre_cut ) continue;
            //float dist =  abs(db_x[in]-q_x) + abs(db_px[in]-q_px) + abs(db_py[in]-q_py) + abs(db_pz[in]-q_pz) ;
            float dist =  (db_x[in]-q_x)*(db_x[in]-q_x) + (db_px[in]-q_px)*(db_px[in]-q_px) + (db_py[in]-q_py)*(db_py[in]-q_py) + (db_pz[in]-q_pz)*(db_pz[in]-q_pz) ;
            if( dist < min_dist ) {min_dist = dist; index = in;}
            //if( min_dist < 1) break;
        }
    }
    //std::cout<<"L2 vector conts_size="<<conts.size()<<",index="<<index<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<",dist="<< min_dist <<std::endl;
    return index; 
}

long FrozenShowerPlugin::Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector<long> & conts) 
{  
    long index=-1;
    float min_dist = 5;
    int conts_size = conts.size();
    for(long k=0 ; k < conts_size ; k++)
    {
        long i = conts.at(k);
        if(abs(db_x[i]-q_x) > min_dist || abs(db_px[i]-q_px) > min_dist || abs(db_py[i]-q_py) > min_dist || abs(db_pz[i]-q_pz) > min_dist ) continue;
        float dist =  abs(db_x[i]-q_x) + abs(db_px[i]-q_px) + abs(db_py[i]-q_py) + abs(db_pz[i]-q_pz) ;
        if( dist < min_dist ) {min_dist = dist; index = i;}
        if( min_dist < 1) break;
    }
    //std::cout<<"abs conts_size="<<conts_size<<",index="<<index<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<",dist="<< min_dist <<std::endl;
    return index; 
}

/*
long FrozenShowerPlugin::Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector<long> & conts) 
{  
    long index=-1;
    float min_dist = 10;
    int conts_size = conts.size();
    for(long k=0 ; k < conts_size ; k++)
    {
        long i = conts.at(k);
        float pre_cut = sqrt(min_dist);
        if(abs(db_x[i]-q_x) > pre_cut || abs(db_px[i]-q_px) > pre_cut || abs(db_py[i]-q_py) > pre_cut || abs(db_pz[i]-q_pz) > pre_cut ) continue;
        //float dist = sqrt( (db_x[i]-q_x)*(db_x[i]-q_x) + (db_px[i]-q_px)*(db_px[i]-q_px) + (db_py[i]-q_py)*(db_py[i]-q_py) + (db_pz[i]-q_pz)*(db_pz[i]-q_pz) );
        float dist =  (db_x[i]-q_x)*(db_x[i]-q_x) + (db_px[i]-q_px)*(db_px[i]-q_px) + (db_py[i]-q_py)*(db_py[i]-q_py) + (db_pz[i]-q_pz)*(db_pz[i]-q_pz) ;
        if( dist < min_dist ) {min_dist = dist; index = i;}
        if( min_dist < 1) break;
    }
    //std::cout<<"conts_size="<<conts_size<<",index="<<index<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<",dist="<< min_dist <<std::endl;
    return index; 
}
*/
/*
long FrozenShowerPlugin::Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const long & tot_evt)
{  
    long index=-1;
    float min_dist = 10;
    for(long i=0 ; i < tot_evt ; i++)
    {
        float pre_cut = sqrt(min_dist);
        if(abs(db_x[i]-q_x) > pre_cut || abs(db_px[i]-q_px) > pre_cut || abs(db_py[i]-q_py) > pre_cut || abs(db_pz[i]-q_pz) > pre_cut ) continue;
        //float dist = sqrt( (db_x[i]-q_x)*(db_x[i]-q_x) + (db_px[i]-q_px)*(db_px[i]-q_px) + (db_py[i]-q_py)*(db_py[i]-q_py) + (db_pz[i]-q_pz)*(db_pz[i]-q_pz) );
        float dist =  (db_x[i]-q_x)*(db_x[i]-q_x) + (db_px[i]-q_px)*(db_px[i]-q_px) + (db_py[i]-q_py)*(db_py[i]-q_py) + (db_pz[i]-q_pz)*(db_pz[i]-q_pz) ;
        if( dist < min_dist ) {min_dist = dist; index = i;}
        //if( min_dist < 1) break;
    }
    //std::cout<<"min_dist="<<min_dist<<",index="<<index<<", q_x="<<q_x<<",q_px="<<q_px<<",q_py="<<q_py<<",q_pz="<<q_pz<<",x="<<db_x[index]<<",px="<<db_px[index]<<",py="<<db_py[index]<<",pz="<<db_pz[index]<<std::endl;
    return index; 
}
*/
faiss::IndexIVFPQ* FrozenShowerPlugin::makeIndex(int dim, const std::string& mc_info, std::vector <float>& database, faiss::IndexFlatL2& coarse_quantizer) {
  std::cout << "Using Frozen shower" << std::endl;

  // dimension of the vectors to index
  // make the index object and train it
  //faiss::IndexFlatL2 coarse_quantizer (dim);
  //std::vector <float> database; 
  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  int N_check = 20000000;
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  float  v_1, v_2, v_3, v_4;
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      v_1 = ::atof(s1.c_str());
      v_2 = ::atof(s2.c_str());
      v_3 = ::atof(s3.c_str());
      v_4 = ::atof(s4.c_str());
      database.push_back(v_1);   
      database.push_back(v_2);   
      database.push_back(v_3);   
      database.push_back(v_4);   
      //if(database.size()/dim > N_check) break;  
  }
  size_t nb = database.size()/dim;
  std::cout<<"nb="<<nb<<std::endl;
  // a reasonable number of centroids to index nb vectors
  int ncentroids = int (4 * sqrt (nb));
  // the coarse quantizer should not be dealloced before the index
  // 4 = nb of bytes per code (d must be a multiple of this)
  // 8 = nb of bits per sub-code (almost always 8)
  faiss::IndexIVFPQ* index = new faiss::IndexIVFPQ(&coarse_quantizer, dim, ncentroids, 4, 8);
  // training
  index->verbose = true;
  index->train (nb, database.data());
  // populating the database_em 
  index->add (nb, database.data());
  std::cout<<"imbalance factor: "<<index->invlists->imbalance_factor()<<std::endl;
  return index;
}


int FrozenShowerPlugin::quireIndex(int nq, int dim, const std::string& mc_info, const faiss::IndexIVFPQ* index) {

  std::vector <float> quire; 
  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  float  v_1, v_2, v_3, v_4;
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      v_1 = ::atof(s1.c_str());
      v_2 = ::atof(s2.c_str());
      v_3 = ::atof(s3.c_str());
      v_4 = ::atof(s4.c_str());
      quire.push_back(v_1);   
      quire.push_back(v_2);   
      quire.push_back(v_3);   
      quire.push_back(v_4);   
      if(quire.size()/dim > nq) break;  
  }
  // searching the database_em
  int k = 2;
  std::vector<faiss::Index::idx_t> nns (k * nq);
  std::vector<float>               dis (k * nq);
  clock_t t0, t1;
  t0 = clock();
  index->search (nq, quire.data(), k, dis.data(), nns.data());
  t1 = clock();
  std::cout << "Hi, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
  for (int i = 0; i < nq; i++) {
      std::cout<< "query "<<i<<": "<<std::endl;
      for (int j = 0; j < k; j++) {
          printf ("%7ld ", nns[j + i * k]);
      }
      printf ("\n     dis: ");
      for (int j = 0; j < k; j++) {
          printf ("%7g ", dis[j + i * k]);
      }
  }
  return 0;
}
