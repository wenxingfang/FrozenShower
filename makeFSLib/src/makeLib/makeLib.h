#ifndef makeLib_h
#define makeLib_h 1


#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"



class makeLib {

    public:
        makeLib(std::string input, std::string output);
        ~makeLib();
        bool configure();               
        bool mutate();    
        bool finish();
        bool isEnd();
        bool clear();
    private:
        long m_total_event;
        long m_processed_event;

        std::vector<double>* m_mc_vertexX;
        std::vector<double>* m_mc_vertexY;
        std::vector<double>* m_mc_vertexZ;
        std::vector<double>* m_mc_Px;
        std::vector<double>* m_mc_Py;
        std::vector<double>* m_mc_Pz;
        std::vector<double>* m_mc_M;
        std::vector<double>* m_mc_Charge; 
        std::vector<double>* m_Hit_x;
        std::vector<double>* m_Hit_y;
        std::vector<double>* m_Hit_z;
        std::vector<double>* m_Hit_E; 
        std::vector<double>* m_HcalHit_x;
        std::vector<double>* m_HcalHit_y;
        std::vector<double>* m_HcalHit_z;
        std::vector<double>* m_HcalHit_E; 

        std::vector< std::vector<double> > m_ECAL_Hit_x; 
        std::vector< std::vector<double> > m_ECAL_Hit_y; 
        std::vector< std::vector<double> > m_ECAL_Hit_z; 
        std::vector< std::vector<double> > m_ECAL_Hit_E; 
        std::vector<int> m_Tag;
        TChain* tree_in;
        TFile* file_out;
        TTree* tree_out;
        TH1D* h_E_n;
        TH1D* h_E_tag0;
        TH1D* h_E_tag1;
        TH1D* h_E_tag2;
        TH1D* h_x_n;
        TH1D* h_x_tag0;
        TH1D* h_x_tag1;
        TH1D* h_x_tag2;
        TH1D* h_theta_n;
        TH1D* h_theta_tag0;
        TH1D* h_theta_tag1;
        TH1D* h_theta_tag2;
        TH2D* h2_x_SumHitE;
        TH2D* h2_x_SumHitE_ori;

};

#endif

