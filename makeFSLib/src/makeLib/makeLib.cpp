#include "makeLib.h"

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
using namespace std;

// add HCALBarrel hit //

makeLib::makeLib(string input, string output){


tree_in = new TChain("evt");
tree_in->Add(input.c_str());
tree_in->SetBranchAddress("m_mc_vertexX", &m_mc_vertexX  );
tree_in->SetBranchAddress("m_mc_vertexY", &m_mc_vertexY  );
tree_in->SetBranchAddress("m_mc_vertexZ", &m_mc_vertexZ  );
tree_in->SetBranchAddress("m_mc_Px", &m_mc_Px       );
tree_in->SetBranchAddress("m_mc_Py", &m_mc_Py       );
tree_in->SetBranchAddress("m_mc_Pz", &m_mc_Pz       );
tree_in->SetBranchAddress("m_mc_M" , &m_mc_M        );
tree_in->SetBranchAddress("m_Hit_x", &m_Hit_x       );
tree_in->SetBranchAddress("m_Hit_y", &m_Hit_y       );
tree_in->SetBranchAddress("m_Hit_z", &m_Hit_z       );
tree_in->SetBranchAddress("m_Hit_E", &m_Hit_E       ); 
tree_in->SetBranchAddress("m_HcalHit_x", &m_HcalHit_x       );
tree_in->SetBranchAddress("m_HcalHit_y", &m_HcalHit_y       );
tree_in->SetBranchAddress("m_HcalHit_z", &m_HcalHit_z       );
tree_in->SetBranchAddress("m_HcalHit_E", &m_HcalHit_E       ); 


file_out = new TFile(output.c_str(),"RECREATE"); 
tree_out = new TTree("evt","tree");
tree_out->Branch("m_ECAL_Hit_x", &m_ECAL_Hit_x       );
tree_out->Branch("m_ECAL_Hit_y", &m_ECAL_Hit_y       );
tree_out->Branch("m_ECAL_Hit_z", &m_ECAL_Hit_z       );
tree_out->Branch("m_ECAL_Hit_E", &m_ECAL_Hit_E       ); 
tree_out->Branch("m_Tag", &m_Tag   ); 

h_E_n       =new TH1D("E_nFS" ,"",500,0,500);          
h_E_tag0    =new TH1D("E_tag0","",500,0,500);
h_E_tag1    =new TH1D("E_tag1","",500,0,500);
h_E_tag2    =new TH1D("E_tag2","",500,0,500);
h_x_n       =new TH1D("x_nFS" ,"",160,0,160);
h_x_tag0    =new TH1D("x_tag0","",160,0,160);
h_x_tag1    =new TH1D("x_tag1","",160,0,160);
h_x_tag2    =new TH1D("x_tag2","",160,0,160);
h_theta_n   =new TH1D("theta_nFS"    ,"",180,0,180);
h_theta_tag0=new TH1D("theta_tag0","",180,0,180);
h_theta_tag1=new TH1D("theta_tag1","",180,0,180);
h_theta_tag2=new TH1D("theta_tag2","",180,0,180);
h2_x_SumHitE=new TH2D("x_SumHitE" ,"",160,1850,2010, 500,0,500);          
h2_x_SumHitE_ori = new TH2D("x_SumHitE_ori" ,"",160,1850,2010, 500,0,500);          

std::cout<<"pi="<< PI << std::endl;

}

makeLib::~makeLib(){
}

bool makeLib::mutate(){
    const int x_bin=160;
    const int theta_bin=180;
    const int E_bin=500;
    const int bins = x_bin*theta_bin*E_bin;// 1850-2009, 0-180, 0-500MeV
    vector< vector< vector<double> > > tmp_Hit_x(bins); 
    vector< vector< vector<double> > > tmp_Hit_y(bins); 
    vector< vector< vector<double> > > tmp_Hit_z(bins); 
    vector< vector< vector<double> > > tmp_Hit_E(bins); 
    vector< vector< int > >            tmp_Tag  (bins); 

    std::cout<<"tmp_Hit_x size:"<<tmp_Hit_x.size()<<",tmp_Tag size="<<tmp_Tag.size()<<std::endl;
    std::cout<<"total events:"<<tree_in->GetEntries()<<std::endl;
    for(int i=0;i<tree_in->GetEntries();i++){
        if(i%1000000==0)std::cout<<"processed:"<<100*i/tree_in->GetEntries()<<"%"<<std::endl;
        tree_in->GetEntry(i);
        if (m_mc_Px->size() != 1) continue;
        int bin_x = int(m_mc_vertexX->at(0)-1850);        
        float mom = sqrt(m_mc_Px->at(0)*m_mc_Px->at(0) + m_mc_Py->at(0)*m_mc_Py->at(0) + m_mc_Pz->at(0)*m_mc_Pz->at(0) );
        float cos_theta = m_mc_Px->at(0)/mom;
        int bin_theta = int(acos(cos_theta)*180/PI);
        int bin_E = int(mom*1000);
        if(bin_E>=500) continue;
        int index = bin_x*180*500 + bin_theta*500 + bin_E;
        if(index<0 || index>bins) continue;
        int tag = 1;
        if(m_Hit_x->size() == 0) tag = 0;
        else if(m_HcalHit_x->size() > 0) tag = 2;
        tmp_Tag.at(index).push_back(tag);
        
        vector<double> tmp_x;
        vector<double> tmp_y;
        vector<double> tmp_z;
        vector<double> tmp_E;
        double tmp_sum_E = 0;
        for(int j=0;j<m_Hit_x->size(); j++)
        {
            tmp_x.push_back(m_Hit_x->at(j));
            tmp_y.push_back(m_Hit_y->at(j));
            tmp_z.push_back(m_Hit_z->at(j));
            tmp_E.push_back(m_Hit_E->at(j));
            tmp_sum_E += m_Hit_E->at(j);
        }
        tmp_Hit_x.at(index).push_back(tmp_x);
        tmp_Hit_y.at(index).push_back(tmp_y);
        tmp_Hit_z.at(index).push_back(tmp_z);
        tmp_Hit_E.at(index).push_back(tmp_E);
        if(bin_E<100) continue; 
        h2_x_SumHitE_ori->Fill(m_mc_vertexX->at(0), tmp_sum_E*1000);//MeV
        /*
        tmp_Hit_x.at(index).push_back(m_Hit_x);
        tmp_Hit_y.at(index).push_back(m_Hit_y);
        tmp_Hit_z.at(index).push_back(m_Hit_z);
        tmp_Hit_E.at(index).push_back(m_Hit_E);
        */
    }
    for(int k=0; k<bins; k++)
    {
        
        if(k%100000==0)std::cout<<"filled:"<<100*k/bins<<"%"<<std::endl;
        clear();
        m_ECAL_Hit_x = tmp_Hit_x.at(k);
        m_ECAL_Hit_y = tmp_Hit_y.at(k);
        m_ECAL_Hit_z = tmp_Hit_z.at(k);
        m_ECAL_Hit_E = tmp_Hit_E.at(k);
        m_Tag        = tmp_Tag  .at(k);
        tree_out->Fill();
        int tag0 = 0;
        int tag1 = 0;
        int tag2 = 0;
        int bin_x = int(k/(E_bin*theta_bin));
        int bin_theta = int((k-bin_x*E_bin*theta_bin)/E_bin);
        int bin_E     = k-bin_x*E_bin*theta_bin-bin_theta*E_bin;
        for(int kk=0;kk<m_Tag.size();kk++)
        {
            if(m_Tag.at(kk)==0) tag0++;
            else if(m_Tag.at(kk)==1) tag1++;
            else if(m_Tag.at(kk)==2) tag2++;
        }
        h_E_n   ->Fill(bin_E, m_Tag.size());
        h_E_tag0->Fill(bin_E, tag0);
        h_E_tag1->Fill(bin_E, tag1);
        h_E_tag2->Fill(bin_E, tag2);
        h_x_n   ->Fill(bin_x, m_Tag.size());
        h_x_tag0->Fill(bin_x, tag0);
        h_x_tag1->Fill(bin_x, tag1);
        h_x_tag2->Fill(bin_x, tag2);
        h_theta_n   ->Fill(bin_theta, m_Tag.size());
        h_theta_tag0->Fill(bin_theta, tag0);
        h_theta_tag1->Fill(bin_theta, tag1);
        h_theta_tag2->Fill(bin_theta, tag2);
        /*
        if(bin_E<100) continue;
        for(int kk=0;kk<m_ECAL_Hit_E.size();kk++)
        {
            double tmp_sumE = 0 ;
            for(int jj=0;jj<m_ECAL_Hit_E.at(kk).size();jj++) tmp_sumE += m_ECAL_Hit_E.at(kk).at(jj);
            h2_x_SumHitE->Fill((1850+bin_x), tmp_sumE*1000);//MeV
        }
        */
    }
    std::cout << "Done for mutate"<< endl;
    return true;
}


bool makeLib::clear(){

    vector< vector<double> >().swap(m_ECAL_Hit_x     );
    vector< vector<double> >().swap(m_ECAL_Hit_y     );
    vector< vector<double> >().swap(m_ECAL_Hit_z     );
    vector< vector<double> >().swap(m_ECAL_Hit_E     ); 
    vector<int>().swap(m_Tag   );


    return true;
}

bool makeLib::isEnd(){
return false;
}

bool makeLib::configure(){
return true;
}

bool makeLib::finish(){
file_out->cd();
tree_out->Write();

h_E_n       ->Write();
h_E_tag0    ->Write();
h_E_tag1    ->Write();
h_E_tag2    ->Write();
h_x_n       ->Write();
h_x_tag0    ->Write();
h_x_tag1    ->Write();
h_x_tag2    ->Write();
h_theta_n   ->Write();
h_theta_tag0->Write();
h_theta_tag1->Write();
h_theta_tag2->Write();
h2_x_SumHitE->Write();
h2_x_SumHitE_ori->Write();


file_out->Close();
return true;
}
