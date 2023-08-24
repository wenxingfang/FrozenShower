#include "makeLib.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "Compression.h"

#include <map>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <cmath> 
#include "my_utils.h"

#ifndef PI
#define PI acos(-1)
#endif
#ifndef DEBUG
#define DEBUG true
#endif
using namespace std;
// add electron and positron
// check phi, theta, E bins, fixed x
// reduce phi bin(-25 to 25, 50 bins) and theta bin (50 to 90, 40 bins)
// increase E bin from 100 MeV to 1 GeV, 900 bin

makeLib::makeLib(string input, string output, string McInfoOut, string s_pdg, string s_phi){

if(s_pdg == "e-") m_pdg = 11;
else if(s_pdg == "e+") m_pdg = -11;
else if(s_pdg == "gamma") m_pdg = 22;

if     (s_phi == "phi+") m_phi = 0;
else if(s_phi == "phi-") m_phi = 1;
else  m_phi = 2;

m_cout = McInfoOut;

tree_in = new TChain("evt");
tree_in->Add(input.c_str());
tree_in->SetBranchAddress("m_mc_vertexX", &m_mc_vertexX  );
tree_in->SetBranchAddress("m_mc_vertexY", &m_mc_vertexY  );
tree_in->SetBranchAddress("m_mc_vertexZ", &m_mc_vertexZ  );
tree_in->SetBranchAddress("m_mc_Px", &m_mc_Px       );
tree_in->SetBranchAddress("m_mc_Py", &m_mc_Py       );
tree_in->SetBranchAddress("m_mc_Pz", &m_mc_Pz       );
tree_in->SetBranchAddress("m_mc_M" , &m_mc_M        );
tree_in->SetBranchAddress("m_mc_Pdg", &m_mc_Pdg     );
tree_in->SetBranchAddress("m_Hit_x", &m_Hit_x       );
tree_in->SetBranchAddress("m_Hit_y", &m_Hit_y       );
tree_in->SetBranchAddress("m_Hit_z", &m_Hit_z       );
tree_in->SetBranchAddress("m_Hit_E", &m_Hit_E       ); 
tree_in->SetBranchAddress("m_HcalHit_x", &m_HcalHit_x       );
tree_in->SetBranchAddress("m_HcalHit_y", &m_HcalHit_y       );
tree_in->SetBranchAddress("m_HcalHit_z", &m_HcalHit_z       );
tree_in->SetBranchAddress("m_HcalHit_E", &m_HcalHit_E       ); 

tree_in->SetBranchAddress("m_ID_S", &m_ID_S      );
tree_in->SetBranchAddress("m_ID_M", &m_ID_M      );
tree_in->SetBranchAddress("m_ID_I", &m_ID_I      );
tree_in->SetBranchAddress("m_ID_J", &m_ID_J      );
tree_in->SetBranchAddress("m_ID_K", &m_ID_K      );

file_out = new TFile(output.c_str(),"RECREATE"); 
//file_out->SetCompressionLevel(ROOT::RCompressionSetting::ELevel::kUncompressed);
file_out->SetCompressionLevel(0);
tree_out = new TTree("evt","tree");
tree_out->Branch("m_ECAL_Hit_x", &m_ECAL_Hit_x       );
tree_out->Branch("m_ECAL_Hit_y", &m_ECAL_Hit_y       );
tree_out->Branch("m_ECAL_Hit_z", &m_ECAL_Hit_z       );
tree_out->Branch("m_ECAL_Hit_E", &m_ECAL_Hit_E       ); 
tree_out->Branch("m_Tag", &m_Tag   ); 
//tree_out->Branch("m_ECAL_ID_S", &m_ECAL_ID_S      );
//tree_out->Branch("m_ECAL_ID_M", &m_ECAL_ID_M      );
tree_out->Branch("m_ECAL_ID_I", &m_ECAL_ID_I      );
//tree_out->Branch("m_ECAL_ID_J", &m_ECAL_ID_J      );
tree_out->Branch("m_ECAL_ID_K", &m_ECAL_ID_K      );

h_mc_mom    =new TH1D("h_mc_mom","",1100,0,1100);          
h_mc_theta  =new TH1D("h_mc_theta","",180,0  ,180);          
h_mc_phi    =new TH1D("h_mc_phi",""  ,180,-90,90);          
h_mc_v_x    =new TH1D("h_mc_v_x","",160,1850,2010);          
h_mc_v_y    =new TH1D("h_mc_v_y","",120,-600,600);          
h_mc_v_z    =new TH1D("h_mc_v_z","",400,-2000,2000);          
h_E_n       =new TH1D("E_nFS" ,"",500,0,500);          
h_E_tag0    =new TH1D("E_tag0","",500,0,500);
h_E_tag1    =new TH1D("E_tag1","",500,0,500);
h_E_tag2    =new TH1D("E_tag2","",500,0,500);
h_phi_n       =new TH1D("phi_nFS" ,"",180,0,180);
h_phi_tag0    =new TH1D("phi_tag0","",180,0,180);
h_phi_tag1    =new TH1D("phi_tag1","",180,0,180);
h_phi_tag2    =new TH1D("phi_tag2","",180,0,180);
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
    std::ofstream mycout(m_cout.c_str());
    std::cout<<"total events:"<<tree_in->GetEntries()<<std::endl;
    for(long long i=0;i<tree_in->GetEntries();i++){
        if(i%1000000==0)std::cout<<"processed:"<<100*i/tree_in->GetEntries()<<"%"<<std::endl;
        tree_in->GetEntry(i);
        if (m_mc_Px->size() != 1) continue;
        if (m_mc_Pdg->at(0) != m_pdg || m_mc_Px->at(0) < 0 || m_mc_Pz->at(0) < 0 ) continue;
        if(m_phi==0 && m_mc_Py->at(0) < 0) continue;
        if(m_phi==1 && m_mc_Py->at(0) > 0) continue;
        mycout << float(m_mc_vertexX->at(0)) <<" " << float(m_mc_Px->at(0)*1000) << " " << float(m_mc_Py->at(0)*1000) <<" " << float(m_mc_Pz->at(0)*1000) <<std::endl;      
        m_Tag = 1;
        if(m_Hit_x->size() == 0) m_Tag = 0;
        else if(m_HcalHit_x->size() > 0) m_Tag = 2;
        clear();
        for(int j=0;j<m_Hit_x->size(); j++)
        {
            m_ECAL_Hit_x.push_back(m_Hit_x->at(j) - m_mc_vertexX->at(0));
            m_ECAL_Hit_y.push_back(m_Hit_y->at(j) - m_mc_vertexY->at(0));
            m_ECAL_Hit_z.push_back(m_Hit_z->at(j) - m_mc_vertexZ->at(0));
            m_ECAL_Hit_E.push_back(m_Hit_E->at(j));
            m_ECAL_ID_S.push_back(m_ID_S->at(j));
            int new_J = m_ID_J->at(j) - int(m_mc_vertexZ->at(0)/10)%90;//
            if (new_J < 0){ new_J = new_J + 90;}
            else if (new_J > 89){ new_J = new_J - 90;}
            m_ECAL_ID_M.push_back(m_ID_M->at(j));//keep unchanged now
            m_ECAL_ID_J.push_back(new_J);//need more check
            m_ECAL_ID_I.push_back(m_ID_I->at(j) + int(m_mc_vertexY->at(0)/10));
            m_ECAL_ID_K.push_back(m_ID_K->at(j));
        }
        tree_out->Fill();
    }
    mycout.close();
    std::cout << "Done for mutate"<< endl;
    return true;
}


bool makeLib::clear(){

    vector<float> ().swap(m_ECAL_Hit_x     );
    vector<float> ().swap(m_ECAL_Hit_y     );
    vector<float> ().swap(m_ECAL_Hit_z     );
    vector<float> ().swap(m_ECAL_Hit_E     ); 
    vector<int> ().swap(m_ECAL_ID_S     ); 
    vector<int> ().swap(m_ECAL_ID_M     ); 
    vector<int> ().swap(m_ECAL_ID_I     ); 
    vector<int> ().swap(m_ECAL_ID_J     ); 
    vector<int> ().swap(m_ECAL_ID_K     ); 

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

h_mc_mom       ->Write();
h_mc_theta     ->Write(); 
h_mc_phi       ->Write(); 
h_mc_v_x       ->Write(); 
h_mc_v_y       ->Write(); 
h_mc_v_z       ->Write(); 

h_E_n       ->Write();
h_E_tag0    ->Write();
h_E_tag1    ->Write();
h_E_tag2    ->Write();
h_phi_n       ->Write();
h_phi_tag0    ->Write();
h_phi_tag1    ->Write();
h_phi_tag2    ->Write();
h_theta_n   ->Write();
h_theta_tag0->Write();
h_theta_tag1->Write();
h_theta_tag2->Write();
h2_x_SumHitE->Write();
h2_x_SumHitE_ori->Write();


file_out->Close();
return true;
}
