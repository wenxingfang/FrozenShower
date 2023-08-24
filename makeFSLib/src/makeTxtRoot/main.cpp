#include "makeLib.h"

#include <iostream>
#include <map>

using namespace std;

int main(int argc, char* argv[]){

for(int i=0;i<argc;i++) cout<<argv[i]<<endl;

std::string generatorName = argv[1];
std::string outputName    = argv[2];
std::string mc_info_output= argv[3];
std::string mc_pdg= argv[4];
std::string s_phi= argv[5];
std::string s_region= argv[6];//barrel or endcap

makeLib* maker  = new makeLib(generatorName, outputName, mc_info_output, mc_pdg, s_phi, s_region);

maker->mutate();
maker->finish();
cout<<"done"<<endl;
return 0;
}
