#include "checkLib.h"

#include <iostream>
#include <map>

using namespace std;

int main(int argc, char* argv[]){

for(int i=0;i<argc;i++) cout<<argv[i]<<endl;

std::string outputName    = argv[1];

makeLib* maker  = new makeLib(outputName);

maker->mutate();
maker->finish();
cout<<"done"<<endl;
return 0;
}
