#include "SLCIORdr.h"

#include <iostream>
#include <map>

using namespace std;

int main(int argc, char* argv[]){

for(int i=0;i<argc;i++) cout<<argv[i]<<endl;

std::string generatorName = argv[1];
std::string outputName    = argv[2];
std::string str_is_gun    = argv[3];
//std::string particle      = argv[4];

bool is_gun = false;
if(str_is_gun=="gun")
{
    cout<<"is gun sample"<<endl;
    is_gun = true;
}
/*
map<string, int> pidMap ;
pidMap.insert(map<string, int>::value_type("e",11));
pidMap.insert(map<string, int>::value_type("gamma",22));
pidMap.insert(map<string, int>::value_type("pi-",-211));
map<string, int>::iterator it;
it = pidMap.find(particle); 
if( it== pidMap.end())
{
    cout<<"wrong particle name"<<endl;
    return 0;
}
std::cout<<"pid="<<it->second<<std::endl;
*/

//SLCIORdr* reader  = new SLCIORdr(generatorName, outputName, is_gun, it->second);
SLCIORdr* reader  = new SLCIORdr(generatorName, outputName, is_gun, 0);

int N=0;
while(reader->mutate() && N<1e6){N++;}
//while(reader->mutate() && N<10){N++;}
reader->finish();
cout<<"done"<<endl;
return 0;
}
