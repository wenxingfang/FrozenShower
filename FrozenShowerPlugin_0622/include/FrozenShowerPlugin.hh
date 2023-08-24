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

#include "VEncoder.hh"
#include "Plugin.hh"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <map>

#include <faiss/IndexIVFFlat.h>
#include <faiss/IndexIVFPQ.h>
#include <faiss/IndexFlat.h>
#include <faiss/index_io.h>

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
  void getID_x_y_z(const int& S, const int& M, const int& I, const int& J, const int& K, int& id, float& x, float& y, float& z);
  int split(vector<string>& res, const string& str, const string& delim) ;
  int quireIndex(int nq, int dim, const std::string& mc_info, const faiss::IndexIVFPQ* index) ;
  faiss::IndexIVFPQ* makeIndex(int dim, const std::string& mc_info, std::vector <float>& database, faiss::IndexFlatL2& coarse_quantizer) ;
  long Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector<long> & conts);
  long Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const long & tot_evt);
  long Search(const float& q_x, const float& q_px, const float& q_py, const float& q_pz, const float* db_x, const float* db_px, const float* db_py , const float* db_pz, const std::vector< std::vector<long> > & conts);

protected:
  FrozenShowerPlugin();

private:
  TChain* tree_in;
  map<int, TChain*> x_chain_map;
  TFile* file_out;
  TTree* tree_out;
  double x_min ;
  double x_max ;//mm
  double theta_min ;
  double theta_max ;
  double phi_min ;
  double phi_max ;
  double Mom_min ;//MeV
  double Mom_max ;//MeV
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
  int m_pass0;
  int m_pass1;
  int m_pass2;
  float m_edep_sum;
  //faiss::IndexIVFPQ* m_index_em;
  //faiss::IndexIVFPQ* m_index_ep;
  faiss::IndexIVFFlat* m_index_em;
  faiss::IndexIVFFlat* m_index_ep;
  std::vector <float> m_database_em; 
  std::vector <float> m_database_ep; 
  faiss::IndexFlatL2 coarse_quantizer_em ;
  faiss::IndexFlatL2 coarse_quantizer_ep ;

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
  float* m_em_x ;
  float* m_em_px;
  float* m_em_py;
  float* m_em_pz;
  float* m_ep_x ;
  float* m_ep_px;
  float* m_ep_py;
  float* m_ep_pz;
  std::vector<long> m_em_cont_x1900L_pyP;
  std::vector<long> m_em_cont_x1900L_pyM;
  std::vector<long> m_em_cont_x1900H_pyP;
  std::vector<long> m_em_cont_x1900H_pyM;
  std::vector<long> m_ep_cont_x1900L_pyP;
  std::vector<long> m_ep_cont_x1900L_pyM;
  std::vector<long> m_ep_cont_x1900H_pyP;
  std::vector<long> m_ep_cont_x1900H_pyM;

  std::vector< std::vector<long> >  m_em_x_vec;
  std::vector< std::vector<long> >  m_ep_x_vec;
};  

#endif // FrozenShowerPlugin_hh
