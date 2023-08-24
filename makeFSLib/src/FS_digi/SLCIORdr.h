#ifndef SLCIORdr_h
#define SLCIORdr_h 1


#include "lcio.h"
#include "LCIOSTLTypes.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCIO.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "IO/LCReader.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "UTIL/CellIDDecoder.h"


#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"


using namespace lcio;

class SLCIORdr {

    public:
        SLCIORdr(std::string name, std::string output, bool is_gun, int pid);
        ~SLCIORdr();
        bool configure();               
        bool mutate();    
        bool getMCinfo(LCCollection* lcCol, std::vector<int>& mc_Pdg, std::vector<double>& mc_Charge, std::vector<int>& mc_genStatus, std::vector<int>& mc_simStatus, std::vector<double>& mc_vertexX, std::vector<double>& mc_vertexY, std::vector<double>& mc_vertexZ, std::vector<double>& mc_Px, std::vector<double>& mc_Py, std::vector<double>& mc_Pz, std::vector<double>& mc_M , std::vector<int>& mc_np, std::vector<int>& mc_nd, std::vector<double>& mc_pHitx, std::vector<double>& mc_pHity, std::vector<double>& mc_pHitz, std::vector<double>& mc_pHit_theta, std::vector<double>& mc_pHit_phi, std::vector<double>& mc_pHit_rotated, LCCollection* SimHitCol, CellIDDecoder<SimCalorimeterHit> & SimHit_idDecoder, LCCollection* DigiHitCol, CellIDDecoder<CalorimeterHit> & DigiHit_idDecoder, std::vector< std::vector<double> >& Hits, std::vector<double>& mc_pHit_dz, std::vector<double>& mc_pHit_dy, LCCollection* DigiHcalCol, std::vector< std::vector<double> >& HcalHits);
        void getMCinfo0(LCCollection* lcCol, std::vector<int>& mc_Pdg, std::vector<double>& mc_Charge, std::vector<int>& mc_genStatus, std::vector<int>& mc_simStatus, std::vector<double>& mc_vertexX, std::vector<double>& mc_vertexY, std::vector<double>& mc_vertexZ, std::vector<double>& mc_Px, std::vector<double>& mc_Py, std::vector<double>& mc_Pz, std::vector<double>& mc_M);
        bool getDigiHitInfo(CellIDDecoder<CalorimeterHit> & Hit_idDecoder, LCCollection* Col, std::vector<int>& vec_ID0, std::vector<int>& vec_ID1, std::vector<double>& vec_Hit_x, std::vector<double>& vec_Hit_y, std::vector<double>& vec_Hit_z, std::vector<double>& vec_Hit_E , double& Hit_cm_x, double& Hit_cm_y, double& Hit_cm_z, double& Hit_tot_e, std::vector<int>& vec_ID_S, std::vector<int>& vec_ID_M, std::vector<int>& vec_ID_I, std::vector<int>& vec_ID_J, std::vector<int>& vec_ID_K);
        bool getSimHitInfo(CellIDDecoder<CalorimeterHit> & Hit_idDecoder, LCCollection* Col, std::vector<int>& vec_ID0, std::vector<int>& vec_ID1, std::vector<double>& vec_Hit_x, std::vector<double>& vec_Hit_y, std::vector<double>& vec_Hit_z, std::vector<double>& vec_Hit_E , double& Hit_cm_x, double& Hit_cm_y, double& Hit_cm_z, double& Hit_tot_e, std::vector<int>& vec_ID_S, std::vector<int>& vec_ID_M, std::vector<int>& vec_ID_I, std::vector<int>& vec_ID_J, std::vector<int>& vec_ID_K);
        bool getSimHcalHitInfo(LCCollection* Col, std::vector<double>& vec_Hit_x, std::vector<double>& vec_Hit_y, std::vector<double>& vec_Hit_z, std::vector<double>& vec_Hit_E);
        bool getHits(LCCollection* SimCol, CellIDDecoder<SimCalorimeterHit>& idDecoderSim, LCCollection* DigiCol, CellIDDecoder<CalorimeterHit>& idDecoderDigi,float x, float y, float z, std::vector<double>& Hits, float& cell_x, float& cell_y, float& cell_z);
        bool getHcalHits(LCCollection* DigiCol, const float& x, const float& y, const float& z, std::vector<double>& Hits, const float& r_min, const float& r_max, const int& r_bin, const float& z_half, const int& z_bin, const float& phi_half, const int& phi_bin);
        bool finish();
        bool isEnd();
        bool clear();
    private:
        IO::LCReader* m_slcio_rdr;
        std::string sim_digi;
        long m_total_event;
        long m_processed_event;
        double m_HitCm_x   ;
        double m_HitCm_y   ;
        double m_HitCm_z   ;
        double m_HitEn_tot ;

        std::vector<double> m_HitFirst_x ;
        std::vector<double> m_HitFirst_y ;
        std::vector<double> m_HitFirst_z ;
        std::vector<double> m_HitFirst_vtheta ;
        std::vector<double> m_HitFirst_vphi   ;
        std::vector<double> m_phi_rotated     ;

        std::vector<double> m_mc_pHitx        ;
        std::vector<double> m_mc_pHity        ;
        std::vector<double> m_mc_pHitz        ;
        std::vector<double> m_mc_pHit_theta   ;
        std::vector<double> m_mc_pHit_phi     ;
        std::vector<double> m_mc_pHit_rotated ;
        std::vector<double> m_mc_pHit_dz ;
        std::vector<double> m_mc_pHit_dy ;


        std::vector<double> m_mc_vertexX;
        std::vector<double> m_mc_vertexY;
        std::vector<double> m_mc_vertexZ;
        std::vector<double> m_mc_Px;
        std::vector<double> m_mc_Py;
        std::vector<double> m_mc_Pz;
        std::vector<double> m_mc_M;
        std::vector<double> m_mc_Charge; 
        std::vector<int   > m_Hit_region;
        std::vector<double> m_Hit_x;
        std::vector<double> m_Hit_y;
        std::vector<double> m_Hit_z;
        std::vector<double> m_Hit_E; 
        std::vector<double> m_HcalHit_x;
        std::vector<double> m_HcalHit_y;
        std::vector<double> m_HcalHit_z;
        std::vector<double> m_HcalHit_E; 
        std::vector< std::vector<double> > m_Hits; 
        std::vector< std::vector<double> > m_HcalHits; 
        std::vector<int>   m_mc_Pdg, m_mc_genStatus, m_mc_simStatus, m_mc_np, m_mc_nd;
        std::vector<int>  m_Hit_ID0,m_Hit_ID1, m_ID_S, m_ID_M, m_ID_I, m_ID_J, m_ID_K;
        TFile* file_out;
        TTree* tree_out;
        int m_pid;
        bool m_is_gun;

};

#endif

