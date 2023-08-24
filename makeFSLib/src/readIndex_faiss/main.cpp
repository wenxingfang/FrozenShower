#include <iostream>
#include <map>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <fstream>


#include <faiss/IndexHNSW.h>
#include <faiss/IndexIVFPQ.h>
#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/index_io.h>


int quireIndex(int nq, int dim, const std::string& mc_info, faiss::IndexIVFPQ* index) ;

using namespace std;
using namespace faiss;

template<typename T>
int quireIndex(int nq, int dim, const std::string& mc_info, T index) ;

int main(int argc, char* argv[]){

for(int i=0;i<argc;i++) cout<<argv[i]<<endl;

std::string inputName     = argv[1];
std::string inputCheckName     = argv[2];
std::cout<<"input index="<<inputName<<std::endl;
std::cout<<"input check file="<<inputCheckName<<std::endl;
//faiss::IndexIVFPQ* index = dynamic_cast<faiss::IndexIVFPQ*> (read_index(inputName.c_str()));
faiss::IndexIVFFlat* index = dynamic_cast<faiss::IndexIVFFlat*> (read_index(inputName.c_str()));
index->nprobe = 1;//default
//faiss::IndexHNSW* index = dynamic_cast<faiss::IndexHNSW*> (read_index(inputName.c_str()));

int s =  quireIndex(1000, 4, inputCheckName, index);
cout<<"done"<<endl;
return 0;
}


template<typename T>
//int quireIndex(int nq, int dim, const std::string& mc_info, faiss::IndexIVFPQ* index) {
int quireIndex(int nq, int dim, const std::string& mc_info, T index) {

  std::vector <float> quire; 
  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  float  v_1, v_2, v_3, v_4;
  float bias = 0;//0.1;
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      v_1 = ::atof(s1.c_str()) + bias;
      v_2 = ::atof(s2.c_str()) + bias;
      v_3 = ::atof(s3.c_str()) + bias;
      v_4 = ::atof(s4.c_str()) + bias;
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
  //clock_t t0, t1;
  //t0 = clock();
  //index->nprobe = 10;
  index->search (nq, quire.data(), k, dis.data(), nns.data());
  //t1 = clock();
  //std::cout << "Hi, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
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

  clock_t t0, t1;
  for(int i=0; i< nq; i++)
  {
      std::cout<< "query each "<<i<<": "<<std::endl;
      std::vector<float> tmp;
      tmp.push_back(quire.at(i*4  )); 
      tmp.push_back(quire.at(i*4+1)); 
      tmp.push_back(quire.at(i*4+2)); 
      tmp.push_back(quire.at(i*4+3)); 
      std::vector<faiss::Index::idx_t> tmp_nns (k * 1);
      std::vector<float>               tmp_dis (k * 1);
      t0 = clock();
      index->search (1, tmp.data(), k, tmp_dis.data(), tmp_nns.data());
      t1 = clock();
      std::cout << "Hi, dt="<<(t1-t0) / (double) CLOCKS_PER_SEC<<",t0="<<t0<<",t1="<<t1<<std::endl;
      for (int j = 0; j < k; j++) {
          printf ("id=%7ld ", tmp_nns[j]);
          printf ("dis=%7g ", tmp_dis[j]);
      }
  }

  return 0;
}
