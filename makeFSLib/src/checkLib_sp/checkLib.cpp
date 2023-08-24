#include "checkLib.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "Compression.h"

#include <map>
#include <sstream>
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


namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

makeLib::makeLib(string output){

  file_out = new TFile(output.c_str(),"RECREATE");


  int x_array[] = {
1870,
1871,
1872,
1873,
1874,
1875,
1876,
1877,
1878,
1879,
1880,
1881,
1882,
1883,
1884,
1885,
1886,
1887,
1888,
1889,
1890,
1891,
1892,
1893,
1894,
1895,
1896,
1897,
1898,
1899,
1900,
1901,
1902,
1903,
1904,
1905,
1906,
1907,
1908,
1909,
1910,
1911,
1912,
1913,
1914,
1915,
1916,
1917,
1918,
1919,
1920,
1921,
1922,
1923,
1924,
1925,
1926,
1927,
1928,
1929,
1930,
1931,
1932,
1933,
1934,
1935,
1936,
1937,
1938,
1939,
1940,
1941,
1942,
1943,
1944,
1945,
1946,
1947,
1948,
1949,
1950,
1951,
1952,
1953,
1954,
1955,
1956,
1957,
1958,
1959,
1960,
1961,
1962,
1963,
1964,
1965,
1966,
1967,
1968,
1969,
1970,
1971,
1972,
1973,
1974,
1975,
1976,
1977,
1978,
1979,
1980,
1981,
1982,
1983,
1984,
1985,
1986,
1987,
1988,
1989,
1990,
1991,
1992,
1993,
1994,
1995,
1996,
1997
};
  for(unsigned int i =0 ; i < (sizeof(x_array) / sizeof(x_array[0])); i++)
  {
     TChain* tmp_tree = new TChain("evt");
     std::string x_s = patch::to_string(x_array[i]);
     std::string file = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0510/lib_x"+x_s+".root";
     tmp_tree->Add(file.c_str());
     tmp_tree->SetBranchAddress("m_ECAL_Hit_x", &m_ECAL_Hit_x       );
     tmp_tree->SetBranchAddress("m_ECAL_Hit_y", &m_ECAL_Hit_y       );
     tmp_tree->SetBranchAddress("m_ECAL_Hit_z", &m_ECAL_Hit_z       );
     tmp_tree->SetBranchAddress("m_ECAL_Hit_E", &m_ECAL_Hit_E       );
     tmp_tree->SetBranchAddress("m_Tag"       , &m_Tag              );
     tmp_tree->SetBranchAddress("m_ECAL_ID_I" , &m_ECAL_ID_I        );
     tmp_tree->SetBranchAddress("m_ECAL_ID_K" , &m_ECAL_ID_K        );
     x_chain_map.insert(pair<int,  TChain*>(x_array[i] ,tmp_tree));    
     total_entry = tmp_tree->GetEntries();
  }


h_ratio_empty_bin    =new TH1D("h_ratio_empty_bin"  ,"",160,1850,2010);          
h_ratio_zero_vhit     =new TH1D("h_ratio_zero_vhit" ,"",160,1850,2010);          

std::cout<<"pi="<< PI << std::endl;

}

makeLib::~makeLib(){
}

bool makeLib::mutate(){
    //const int pdg_bin=2;
    //const int phi_bin=50;
    //const int theta_bin=40;
    //const int E_bin=900;
    //const int bins = pdg_bin*phi_bin*theta_bin*E_bin;
    for(std::map<int, TChain*>::iterator iter=x_chain_map.begin(); iter != x_chain_map.end(); iter++)
    {
        int empty_bin = 0;
        int all_vhit = 0;
        int zero_vhit = 0;
        int x_tmp = iter->first;
        int entries = iter->second->GetEntries();
        for(int i=0; i< entries; i++)
        {
            iter->second->GetEntry(i);
            if(m_Tag->size() == 0) empty_bin ++; 
            else
            {
                for(int j=0; j<m_Tag->size(); j++)
                {   all_vhit ++;
                    if(m_Tag->at(j) == 0) zero_vhit ++;
                }
            }
        }
        std::cout<<"x_tmp="<<x_tmp<<",tot="<<entries<<",empty_bin="<<empty_bin<<",all_vhit="<<all_vhit<<",zero_vhit="<<zero_vhit<<std::endl;
        float ratio_empty_bin = entries  != 0 ? float(empty_bin)/entries : 0 ;
        float ratio_zero_vhit = all_vhit != 0 ? float(zero_vhit)/all_vhit : 0 ;
        h_ratio_empty_bin->Fill(x_tmp,ratio_empty_bin) ;
        h_ratio_zero_vhit->Fill(x_tmp,ratio_zero_vhit) ;
    }

    std::cout << "Done for mutate"<< endl;
    return true;
}


bool makeLib::clear(){


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

h_ratio_empty_bin    ->Write();
h_ratio_zero_vhit    ->Write(); 

file_out->Close();
return true;
}
