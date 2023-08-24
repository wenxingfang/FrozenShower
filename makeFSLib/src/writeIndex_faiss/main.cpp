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

#include <omp.h>

#include <faiss/IndexHNSW.h>
#include <faiss/index_factory.h>
#include <faiss/IndexIVFPQ.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/IndexFlat.h>
#include <faiss/index_io.h>

faiss::IndexIVFPQ* makeIndex(int dim, const std::string& mc_info, std::vector <float>& database, faiss::IndexFlatL2& coarse_quantizer) ;
faiss::IndexIVFFlat* makeIndex(int dim, const std::string& mc_info, std::vector <float>& database) ;
faiss::IndexHNSW* makeIndexHNSW(int dim, const std::string& mc_info, std::vector <float>& database) ;

using namespace std;

int main(int argc, char* argv[]){

for(int i=0;i<argc;i++) cout<<argv[i]<<endl;

std::string inputName     = argv[1];
std::string outputName    = argv[2];

omp_set_num_threads(1);

  int dim = 4;
  //int dim_ep = 4;
  faiss::IndexFlatL2 coarse_quantizer (dim);
  //faiss::IndexFlatL2 coarse_quantizer_ep (dim_ep);
  //std::string mc_info_em = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0514//lib_e-.txt";
  //std::string mc_info_ep = "/cefs/higgs/wxfang/cepc/FS/e_startpoint_0514//lib_e+.txt";
  //m_index_em = makeIndex(dim_em, mc_info_em, m_database_em, coarse_quantizer_em) ;
  std::vector <float> m_database;
  //faiss::IndexIVFPQ* m_index = makeIndex(dim, inputName, m_database, coarse_quantizer) ;
  //faiss::IndexIVFFlat* m_index = makeIndex(dim, inputName, m_database) ;
  faiss::IndexHNSW* m_index = makeIndexHNSW(dim, inputName, m_database) ;
  write_index(m_index, outputName.c_str());


cout<<"done"<<endl;
return 0;
}


faiss::IndexIVFPQ* makeIndex(int dim, const std::string& mc_info, std::vector <float>& database, faiss::IndexFlatL2& coarse_quantizer) {
  std::cout << "Using Frozen shower" << std::endl;

  // dimension of the vectors to index
  // make the index object and train it
  //faiss::IndexFlatL2 coarse_quantizer (dim);
  //std::vector <float> database; 
  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  //int N_check = 100000;
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  float  v_1, v_2, v_3, v_4;
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      v_1 = ::atof(s1.c_str());
      v_2 = ::atof(s2.c_str());
      v_3 = ::atof(s3.c_str());
      v_4 = ::atof(s4.c_str());
      database.push_back(v_1);   
      database.push_back(v_2);   
      database.push_back(v_3);   
      database.push_back(v_4);   
      //if(database.size()/dim > N_check) break;  
  }
  size_t nb = database.size()/dim;
  std::cout<<"nb="<<nb<<std::endl;
  // a reasonable number of centroids to index nb vectors
  int ncentroids = int (4 * sqrt (nb));
  // the coarse quantizer should not be dealloced before the index
  // 4 = nb of bytes per code (d must be a multiple of this)
  // 8 = nb of bits per sub-code (almost always 8)
  faiss::IndexIVFPQ* index = new faiss::IndexIVFPQ(&coarse_quantizer, dim, ncentroids, 4, 8);
  // training
  index->verbose = true;
  index->train (nb, database.data());
  // populating the database_em 
  index->add (nb, database.data());
  std::cout<<"imbalance factor: "<<index->invlists->imbalance_factor()<<std::endl;
  return index;
}


faiss::IndexIVFFlat* makeIndex(int dim, const std::string& mc_info, std::vector <float>& database) {
  
  faiss::IndexIVFFlat* index = dynamic_cast<faiss::IndexIVFFlat*> (faiss::index_factory(dim, "IVF65536_HNSW32,Flat"));
  //faiss::IndexHNSWFlat* index = new faiss::IndexHNSWFlat(dim,32);
  // dimension of the vectors to index
  // make the index object and train it
  //faiss::IndexFlatL2 coarse_quantizer (dim);
  //std::vector <float> database; 
  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  int N_check = 5000000;
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  float  v_1, v_2, v_3, v_4;
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      v_1 = ::atof(s1.c_str());
      v_2 = ::atof(s2.c_str());
      v_3 = ::atof(s3.c_str());
      v_4 = ::atof(s4.c_str());
      database.push_back(v_1);   
      database.push_back(v_2);   
      database.push_back(v_3);   
      database.push_back(v_4);   
      //if(database.size()/dim > N_check) break;  
  }
  size_t nb = database.size()/dim;
  std::cout<<"nb="<<nb<<std::endl;
  // a reasonable number of centroids to index nb vectors
  //int ncentroids = int (4 * sqrt (nb));
  // the coarse quantizer should not be dealloced before the index
  // 4 = nb of bytes per code (d must be a multiple of this)
  // 8 = nb of bits per sub-code (almost always 8)
  //faiss::IndexIVFPQ* index = new faiss::IndexIVFPQ(&coarse_quantizer, dim, ncentroids, 4, 8);
  //faiss::IndexHNSWFlat* index = dynamic_cast<faiss::IndexHNSWFlat*> (faiss::index_factory(dim, "IVF65536_HNSW32,Flat"));
  // training
  index->verbose = true;
  index->train (nb, database.data());
  // populating the database_em 
  index->add (nb, database.data());
  return index;
}


faiss::IndexHNSW* makeIndexHNSW(int dim, const std::string& mc_info, std::vector <float>& database) {
  
  // dimension of the vectors to index
  // make the index object and train it
  //faiss::IndexFlatL2 coarse_quantizer (dim);
  //std::vector <float> database; 
  std::ifstream m_f;
  m_f.open(mc_info.c_str());
  if(m_f.fail()) { std::cout << "error, can't open "<<mc_info<< std::endl; throw ;}
  int N_check = 10000000;
  std::string sline;//每一行
  std::string s1,s2,s3,s4;
  float  v_1, v_2, v_3, v_4;
  while(getline(m_f,sline)){
      std::istringstream sin(sline);
      sin>>s1>>s2>>s3>>s4;
      v_1 = ::atof(s1.c_str());
      v_2 = ::atof(s2.c_str());
      v_3 = ::atof(s3.c_str());
      v_4 = ::atof(s4.c_str());
      database.push_back(v_1);   
      database.push_back(v_2);   
      database.push_back(v_3);   
      database.push_back(v_4);   
      //if(database.size()/dim > N_check) break;  
  }
  size_t nb = database.size()/dim;
  std::cout<<"nb="<<nb<<std::endl;
  faiss::IndexHNSW* index = dynamic_cast<faiss::IndexHNSW*> (faiss::index_factory(dim, "HNSW32"));
  index->verbose = true;
  index->add (nb, database.data());
  return index;
}
