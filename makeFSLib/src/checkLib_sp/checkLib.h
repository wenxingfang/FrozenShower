#ifndef checkLib_h
#define checkLib_h 1


#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include <map>


class makeLib {

    public:
        makeLib(std::string output);
        ~makeLib();
        bool configure();               
        bool mutate();    
        bool finish();
        bool isEnd();
        bool clear();
    private:

        std::map<int, TChain*> x_chain_map;
        TFile* file_out;
        TH1D* h_ratio_empty_bin ;
        TH1D* h_ratio_zero_vhit  ;
        int  total_entry;
        std::vector< std::vector<float> >* m_ECAL_Hit_x;
        std::vector< std::vector<float> >* m_ECAL_Hit_y;
        std::vector< std::vector<float> >* m_ECAL_Hit_z;
        std::vector< std::vector<float> >* m_ECAL_Hit_E;
        std::vector< std::vector<int> >* m_ECAL_ID_K;
        std::vector< std::vector<int> >* m_ECAL_ID_I;
        std::vector<int>* m_Tag;
      
};

#endif

