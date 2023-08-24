#include "makeLib.h"

#include <iostream>
#include <map>

using namespace std;

int main(int argc, char* argv[]){

for(int i=0;i<argc;i++) cout<<argv[i]<<endl;

std::string generatorName = argv[1];
std::string outputName    = argv[2];

makeLib* maker  = new makeLib(generatorName, outputName);

maker->mutate();
maker->finish();
cout<<"done"<<endl;
return 0;
}
