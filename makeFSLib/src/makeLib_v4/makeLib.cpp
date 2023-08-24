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
#include "my_utils.h"

#ifndef PI
#define PI acos(-1)
#endif
#ifndef DEBUG
#define DEBUG true
#endif
using namespace std;

// add HCALBarrel hit //
// check phi, theta, E bins, fixed x

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

h_mc_mom    =new TH1D("h_mc_mom","",500,0,500);          
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
    const int phi_bin=180;
    const int theta_bin=90;
    const int E_bin=400;
    const int bins = phi_bin*theta_bin*E_bin;// mom from 100-500MeV
    vector< vector< vector<float> > > tmp_Hit_x(bins); 
    vector< vector< vector<float> > > tmp_Hit_y(bins); 
    vector< vector< vector<float> > > tmp_Hit_z(bins); 
    vector< vector< vector<float> > > tmp_Hit_E(bins); 
    vector< vector< int > >            tmp_Tag  (bins); 
    vector< vector< vector<int> > >    tmp_ID_S (bins); 
    vector< vector< vector<int> > >    tmp_ID_M (bins); 
    vector< vector< vector<int> > >    tmp_ID_I (bins); 
    vector< vector< vector<int> > >    tmp_ID_J (bins); 
    vector< vector< vector<int> > >    tmp_ID_K (bins); 

    std::cout<<"tmp_Hit_x size:"<<tmp_Hit_x.size()<<",tmp_Tag size="<<tmp_Tag.size()<<std::endl;
    std::cout<<"total events:"<<tree_in->GetEntries()<<std::endl;
    for(int i=0;i<tree_in->GetEntries();i++){
        if(i%1000000==0)std::cout<<"processed:"<<100*i/tree_in->GetEntries()<<"%"<<std::endl;
        tree_in->GetEntry(i);
        if (m_mc_Px->size() != 1) continue;
        if (m_mc_Px->at(0) < 0 || m_mc_Pz->at(0) < 0 ) continue;
        //int bin_x = int(m_mc_vertexX->at(0)-1850);        
        float pt = sqrt(m_mc_Px->at(0)*m_mc_Px->at(0) + m_mc_Py->at(0)*m_mc_Py->at(0) );
        float cos_phi = m_mc_Px->at(0)/pt;
        int bin_phi = (m_mc_Py->at(0) > 0) ? int(acos(cos_phi)*180/PI) : 90 + int(acos(cos_phi)*180/PI);
        if(bin_phi>180 || bin_phi<0) continue;
        float mom = sqrt(m_mc_Px->at(0)*m_mc_Px->at(0) + m_mc_Py->at(0)*m_mc_Py->at(0) + m_mc_Pz->at(0)*m_mc_Pz->at(0) );
        float cos_theta = m_mc_Pz->at(0)/mom;
        int bin_theta = int(acos(cos_theta)*180/PI);
        if(bin_theta>90 || bin_theta<0) continue;
        int bin_E = int(mom*1000)-100;
        h_mc_mom  ->Fill( mom*1000 );
        if(bin_E>400 || bin_E<0) continue;
        int index = bin_phi*90*400 + bin_theta*400 + bin_E;
        if(index<0 || index>=bins) continue;
        int tag = 1;
        if(m_Hit_x->size() == 0) tag = 0;
        else if(m_HcalHit_x->size() > 0) tag = 2;
        tmp_Tag.at(index).push_back(tag);
        
        vector<float> tmp_x;
        vector<float> tmp_y;
        vector<float> tmp_z;
        vector<float> tmp_E;
        vector<int>    tmp_S;
        vector<int>    tmp_M;
        vector<int>    tmp_I;
        vector<int>    tmp_J;
        vector<int>    tmp_K;
        float tmp_sum_E = 0;
        for(int j=0;j<m_Hit_x->size(); j++)
        {
            tmp_x.push_back(m_Hit_x->at(j) - m_mc_vertexX->at(0));
            tmp_y.push_back(m_Hit_y->at(j) - m_mc_vertexY->at(0));
            tmp_z.push_back(m_Hit_z->at(j) - m_mc_vertexZ->at(0));
            tmp_E.push_back(m_Hit_E->at(j));
            tmp_sum_E += m_Hit_E->at(j);
            tmp_S.push_back(m_ID_S->at(j));
            int new_J = m_ID_J->at(j) - int(m_mc_vertexZ->at(0)/10)%90;//
            if (new_J < 0){ new_J = new_J + 90;}
            else if (new_J > 89){ new_J = new_J - 90;}
            //new_M = m_ID_M->at(j) - int(m_mc_vertexZ->at(0)/10);
            tmp_M.push_back(m_ID_M->at(j));//keep unchanged now
            //tmp_I.push_back(m_ID_S->at(j));
            tmp_I.push_back(m_ID_I->at(j) + int(m_mc_vertexY->at(0)/10));
            tmp_J.push_back(new_J);
            tmp_K.push_back(m_ID_K->at(j));
        }
        tmp_Hit_x.at(index).push_back(tmp_x);
        tmp_Hit_y.at(index).push_back(tmp_y);
        tmp_Hit_z.at(index).push_back(tmp_z);
        tmp_Hit_E.at(index).push_back(tmp_E);
        h2_x_SumHitE_ori->Fill(m_mc_vertexX->at(0), tmp_sum_E*1000);//MeV
        tmp_ID_S.at(index).push_back(tmp_S);
        tmp_ID_M.at(index).push_back(tmp_M);
        tmp_ID_I.at(index).push_back(tmp_I);
        tmp_ID_J.at(index).push_back(tmp_J);
        tmp_ID_K.at(index).push_back(tmp_K);
        /*
        tmp_Hit_x.at(index).push_back(m_Hit_x);
        tmp_Hit_y.at(index).push_back(m_Hit_y);
        tmp_Hit_z.at(index).push_back(m_Hit_z);
        tmp_Hit_E.at(index).push_back(m_Hit_E);
        */
    }
    for(int k=0; k<bins; k++)
    {
        
        if(k%1000000==0)std::cout<<"filled:"<<100*k/bins<<"%"<<std::endl;
        clear();
        m_ECAL_Hit_x = tmp_Hit_x.at(k);
        m_ECAL_Hit_y = tmp_Hit_y.at(k);
        m_ECAL_Hit_z = tmp_Hit_z.at(k);
        m_ECAL_Hit_E = tmp_Hit_E.at(k);
        m_Tag        = tmp_Tag  .at(k);
        m_ECAL_ID_S  = tmp_ID_S.at(k);
        m_ECAL_ID_M  = tmp_ID_M.at(k);
        m_ECAL_ID_I  = tmp_ID_I.at(k);
        m_ECAL_ID_J  = tmp_ID_J.at(k);
        m_ECAL_ID_K  = tmp_ID_K.at(k);
        tree_out->Fill();
        /*
        int tag0 = 0;
        int tag1 = 0;
        int tag2 = 0;
        int bin_phi   = int(k/(E_bin*theta_bin));
        int bin_theta = int((k-bin_phi*E_bin*theta_bin)/E_bin);
        int bin_E     = k-bin_phi*E_bin*theta_bin-bin_theta*E_bin;
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
        h_phi_n   ->Fill(bin_phi, m_Tag.size());
        h_phi_tag0->Fill(bin_phi, tag0);
        h_phi_tag1->Fill(bin_phi, tag1);
        h_phi_tag2->Fill(bin_phi, tag2);
        h_theta_n   ->Fill(bin_theta, m_Tag.size());
        h_theta_tag0->Fill(bin_theta, tag0);
        h_theta_tag1->Fill(bin_theta, tag1);
        h_theta_tag2->Fill(bin_theta, tag2);
        */
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

    vector< vector<float> >().swap(m_ECAL_Hit_x     );
    vector< vector<float> >().swap(m_ECAL_Hit_y     );
    vector< vector<float> >().swap(m_ECAL_Hit_z     );
    vector< vector<float> >().swap(m_ECAL_Hit_E     ); 
    vector< vector<int> >().swap(m_ECAL_ID_S     ); 
    vector< vector<int> >().swap(m_ECAL_ID_M     ); 
    vector< vector<int> >().swap(m_ECAL_ID_I     ); 
    vector< vector<int> >().swap(m_ECAL_ID_J     ); 
    vector< vector<int> >().swap(m_ECAL_ID_K     ); 
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

h_mc_mom       ->Write();
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
