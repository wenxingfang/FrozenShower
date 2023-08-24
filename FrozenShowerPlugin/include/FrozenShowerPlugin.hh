// ********************************************************
// *            Mokka    http://mokka.in2p3.fr            *
// *    -- A Detailed Geant4 Simulation for the ILC --    *
// ********************************************************
//
//
// "FrozenShowerPlugin" .
//
//

#ifndef FrozenShowerPlugin_hh
#define FrozenShowerPlugin_hh 1

#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "CalHit.hh"
#include "VSensitiveDetector.hh"
#include "SEcalSD02.hh"
#include "VEncoder.hh"
#include "Plugin.hh"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include <map>
#include <time.h>
//#include "RandomEngine.h"
//#include <RandGauss.h>
//#include <faiss/IndexIVFFlat.h>
//#include <faiss/IndexIVFPQ.h>
//#include <faiss/IndexFlat.h>
//#include <faiss/index_io.h>
using namespace std;
class FrozenShowerPlugin: public Plugin
{
public:
  FrozenShowerPlugin(const std::string &name): Plugin(name) {}
  virtual ~FrozenShowerPlugin() {}

  virtual void Init();
  virtual void Exit();

  virtual void BeginOfRunAction(const G4Run *run);
  virtual void EndOfRunAction(const G4Run *run);

  virtual void BeginOfEventAction(const G4Event *evt);
  virtual void EndOfEventAction(const G4Event *evt);

  virtual void PreUserTrackingAction(const G4Track *trk);
  virtual void PostUserTrackingAction(const G4Track *trk);

  virtual void UserSteppingAction(const G4Step *step);
  int makeHitCollection(const G4Step *step, int & _pid, double & _x, double & _y, double & _z, double & _px, double & _py, double & _pz);
  float getPhi(float x, float y);
  int partition(float phi);
  int partition(float x, float y);
  void line_a_b(float x1, float y1, float x2, float y2, float& a, float& b);
  void getID_x_y_z(const int& S, const int& M, const int& I, const int& J, const int& K, int& id, float& x, float& y, float& z);
  int split(vector<string>& res, const string& str, const string& delim) ;
  long Search(const float& low, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector<long> > > & conts, float& en_scale);
  long Search(const float& low, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector<long> > > & conts, float& en_scale, float& b_x, float& b_px, float& b_py, float& b_pz);
  long Search_v1(const float& low, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector<long> > > & conts, float& en_scale, float& b_x, float& b_px, float& b_py, float& b_pz);
  void processMapping( map<int, TChain*>& chain_map, const std::string& x_s, const std::string& file) ;
  void processTxt(const std::string& mc_info, float* & array_0, float*  & array_1, float* & array_2, float* & array_3, bool isBarrel, std::vector< std::vector< std::vector<long> > > & v_v_vec) ;
  void processTxt(const std::string& mc_info, float* & array_0, float*  & array_1, float* & array_2, float* & array_3, bool isBarrel, std::vector< std::vector< std::vector< std::vector< std::vector<long> > > > > & v_v_vec, long & n_tot) ;
  long Search_v1(const float& low_x, const float& low_theta, const float& low_phi, const float& low_mom, const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector< std::vector< std::vector< std::vector<long> > > > > & conts, float& en_scale, float& b_x, float& b_px, float& b_py, float& b_pz, TH1D* h_res);

  float getE_scale(float x, float y); // from -90 to 90
  void createMomHist(TH1D* h1, float* px, float* py, float* pz, long n_tot) ;
  void printHist(const TH1D* h1 );
  void createResHist(const TH1D* h1, TH1D* h_res, float res );
protected:
  FrozenShowerPlugin();

private:
  TChain* tree_in;
  map<int, TChain*> m_barrel_x_chain_map;
  map<int, TChain*> m_endcap_z_chain_map;
  TFile* file_out;
  TTree* tree_out;
  double m_x_min ;
  double m_x_max ;//mm
  double m_theta_min ;
  double m_theta_max ;
  double m_phi_min ;
  double m_phi_max ;
  double m_Mom_min ;//MeV
  double m_Mom_max ;//MeV
  double m_res ;
  double y_boundry ;//mm
  double z_boundry ;//mm
  double outCut ;
  int  total_entry_em;
  int  total_entry_ep;
  std::vector<float>* m_ECAL_Hit_x;
  std::vector<float>* m_ECAL_Hit_y;
  std::vector<float>* m_ECAL_Hit_z;
  std::vector<float>* m_ECAL_Hit_E;
  std::vector<int>* m_ECAL_ID_K;
  std::vector<int>* m_ECAL_ID_I;
  int m_Tag;
  std::vector<int> m_layer;
  std::vector<int> m_return;
  std::vector<double> m_ratio;
  std::vector<int>    m_s_pid ;
  std::vector<double> m_s_x   ;
  std::vector<double> m_s_y   ;
  std::vector<double> m_s_z   ;
  std::vector<double> m_s_px  ;
  std::vector<double> m_s_py  ;
  std::vector<double> m_s_pz  ;

  map<int, map<int, map<int, map<int, map<int, string> > > > >  SMIJK_map_ID_x_y_z;
  VEncoder *theEncoder;
  float width_M ;
  float width_J ;
  bool doSaveSP ;
  std::vector<double> m_point_x;
  std::vector<double> m_point_y;
  std::vector<double> m_point_z;
  std::vector<double> m_mom_x;
  std::vector<double> m_mom_y;
  std::vector<double> m_mom_z;
  std::vector<int> m_pid;
  std::vector<double> m_noFind_point_x;
  std::vector<double> m_noFind_point_y;
  std::vector<double> m_noFind_point_z;
  std::vector<double> m_noFind_mom_x;
  std::vector<double> m_noFind_mom_y;
  std::vector<double> m_noFind_mom_z;
  std::vector<int>    m_noFind_pid;
  bool saveNotFind;

  std::vector<double> m_Find_point_x ;
  std::vector<double> m_Find_point_y ;
  std::vector<double> m_Find_point_z ;
  std::vector<double> m_Find_mom_x   ;
  std::vector<double> m_Find_mom_y   ;
  std::vector<double> m_Find_mom_z   ;
  std::vector<double> m_Find_sum_hit_E;
  std::vector<double> m_Find_sum_hit_e;
  std::vector<int>    m_Find_pid      ;
  bool saveFind;
  std::vector<double> m_FindLib_x ;
  std::vector<double> m_FindLib_mom_x ;
  std::vector<double> m_FindLib_mom_y ;
  std::vector<double> m_FindLib_mom_z ;


  int m_pass0;
  int m_pass1;
  int m_pass2;
  float m_edep_sum;
  //faiss::IndexIVFPQ* m_index_em;
  //faiss::IndexIVFPQ* m_index_ep;
  //faiss::IndexIVFFlat* m_index_em;
  //faiss::IndexIVFFlat* m_index_ep;
  std::vector <float> m_database_em; 
  std::vector <float> m_database_ep; 
  //faiss::IndexFlatL2 coarse_quantizer_em ;
  //faiss::IndexFlatL2 coarse_quantizer_ep ;

  bool checkTime;
  std::vector< std::vector<float> > m_check_em_hits_x;
  std::vector< std::vector<float> > m_check_em_hits_y;
  std::vector< std::vector<float> > m_check_em_hits_z;
  std::vector< std::vector<float> > m_check_em_hits_E;
  std::vector< std::vector<int> > m_check_em_hits_I;
  std::vector< std::vector<int> > m_check_em_hits_K;
  std::vector< int  >               m_check_em_tag   ;
  std::vector< std::vector<float> > m_check_ep_hits_x;
  std::vector< std::vector<float> > m_check_ep_hits_y;
  std::vector< std::vector<float> > m_check_ep_hits_z;
  std::vector< std::vector<float> > m_check_ep_hits_E;
  std::vector< std::vector<int> > m_check_ep_hits_I;
  std::vector< std::vector<int> > m_check_ep_hits_K;
  std::vector< int  >               m_check_ep_tag   ;
  long m_em_total;
  long m_ep_total;
  float* m_em_x_barrel;
  float* m_em_px_barrel;
  float* m_em_py_barrel;
  float* m_em_pz_barrel;
  float* m_ep_x_barrel;
  float* m_ep_px_barrel;
  float* m_ep_py_barrel;
  float* m_ep_pz_barrel;
  //std::vector< std::vector< std::vector<long> > >  m_em_x_barrel_mom_vec;
  //std::vector< std::vector< std::vector<long> > >  m_ep_x_barrel_mom_vec;
  std::vector< std::vector< std::vector< std::vector< std::vector<long> > > > > m_em_barrel_x_theta_phi_mom_vec; 
  std::vector< std::vector< std::vector< std::vector< std::vector<long> > > > > m_ep_barrel_x_theta_phi_mom_vec; 
  float* m_em_z_endcap;
  float* m_em_px_endcap;
  float* m_em_py_endcap;
  float* m_em_pz_endcap;
  float* m_ep_z_endcap;
  float* m_ep_px_endcap;
  float* m_ep_py_endcap;
  float* m_ep_pz_endcap;
  std::vector< std::vector< std::vector<long> > >  m_em_z_endcap_mom_vec;
  std::vector< std::vector< std::vector<long> > >  m_ep_z_endcap_mom_vec;


  TGraphErrors* m_gr_phibin_Emean;
  bool doEnergyTunning;
  bool apply_Region;
  float abs_endcap_z_min ;//mm
  float abs_endcap_z_max ;//mm
  float abs_endcap_inner_x_min ;//mm
  float abs_endcap_inner_y_min ;//mm
  float abs_endcap_outer_r_max ;//mm

  int Endcap_I_index_max ;
  int Endcap_I_index_min ;
  int Endcap_J_index_max ;
  int Endcap_J_index_min ;
  float Endcap_x_max ;//mm
  float Endcap_y_max ;//mm
  float Endcap_I_index_bin;
  float Endcap_J_index_bin;

  float m_barrel_SD_mag ;
  float PI ;
  clock_t t0, t1;
  HitsCollection *barrel_CalCollection;
  SEcalSD02* theSD_barrel;
  TH1D* m_h_mom;
  TH1D* m_h_res;
  float m_global_scale ;
  //HepRandomEngine* m_engine;
  //RandGauss* m_distribution;
public:
    inline double gaussrand()
    {
        static double V1, V2, S;
        static int phase = 0;
        double X;
    
        if ( phase == 0 ) {
            do {
                double U1 = (double)rand() / RAND_MAX;
                double U2 = (double)rand() / RAND_MAX;
    
                V1 = 2 * U1 - 1;
                V2 = 2 * U2 - 1;
                S = V1 * V1 + V2 * V2;
            } while(S >= 1 || S == 0);
    
            X = V1 * sqrt(-2 * log(S) / S);
        } else
            X = V2 * sqrt(-2 * log(S) / S);
    
        phase = 1 - phase;
    
        return X;
    }


};  

inline float FrozenShowerPlugin::getPhi(float x, float y)
{
    if     (x==0 && y>0) return 90;
    else if(x==0 && y<0) return 270;
    else if(x==0 && y==0) return 0;
    float phi = atan(y/x)*180/acos(-1);
    if                 (x<0) phi = phi + 180;
    else if     (x>0 && y<0) phi = phi + 360;
    return phi;
}

inline int FrozenShowerPlugin::partition(float phi) // phi should < 360 and > 0
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

inline int FrozenShowerPlugin::partition(float x, float y)
{
    float a1, b1, a2, b2 ;
    if (x>1850 && x < 2020){
        line_a_b(1850, 750  , 2020, 600 , a1, b1);
        line_a_b(1850, -1000, 2020, -850, a2, b2);
        if ( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1 ) return 1 ;
    }
    if (y < 1850 && x < 2020){
        line_a_b(750 , 1850,2020, 600, a1, b1);
        line_a_b(1000, 1850,2020, 850, a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1 ) return 2;
    }
    if (y < 2020 && y > 1850){
        line_a_b( -750, 1850 , -600, 2020 , a1, b1);
        line_a_b( 1000, 1850 , 850 , 2020 , a2, b2);
        if( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1 ) return 3;
    }
    if (y < 2020 && x > -1850){
        line_a_b(  -1850, 750  , -600 , 2020 , a1, b1);
        line_a_b(  -1850, 1000 , -850 , 2020 , a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1 ) return 4;
    }
    if (x > -2020 && x < -1850){
        line_a_b( -1850, -750 ,  -2020, -600 , a1, b1);
        line_a_b( -1850, 1000 ,  -2020,  850 , a2, b2);
        if( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1 ) return 5;
    }
    if (x > -2020 && y > -1850){
        line_a_b( -750 , -1850 , -2020, -600 , a1, b1);
        line_a_b( -1000, -1850 , -2020, -850 , a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1) return 6;
    }
    if (y > -2020 && y < -1850){
        line_a_b(  750 , -1850 ,  600, -2020 , a1, b1);
        line_a_b( -1000, -1850 , -850, -2020 , a2, b2);
        if( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1) return 7;
    }
    if (y > -2020 && x < 1850){
        line_a_b(  1850, -750  ,  600, -2020 , a1, b1);
        line_a_b(  1850, -1000 ,  850, -2020 , a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1) return 8;
    }
    return -1;
}

inline void FrozenShowerPlugin::line_a_b(float x1, float y1, float x2, float y2, float& a, float& b)
{
    a = (y2-y1)/(x1*y2-x2*y1);
    b = (x1-x2)/(x1*y2-x2*y1);
}
#endif // FrozenShowerPlugin_hh
