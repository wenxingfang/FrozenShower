#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH2F.h"
#include <iostream>
#include <map>

using namespace std;

void line_a_b(float x1, float y1, float x2, float y2, float& a, float& b);
int partition(float x, float y);

int main(int argc, char* argv[]){

for(int i=0;i<argc;i++) cout<<argv[i]<<endl;

TGraph* gr1=new TGraph();
TGraph* gr2=new TGraph();
TGraph* gr3=new TGraph();
TGraph* gr4=new TGraph();
TGraph* gr5=new TGraph();
TGraph* gr6=new TGraph();
TGraph* gr7=new TGraph();
TGraph* gr8=new TGraph();

long tot = 0;
while(1){
    if (tot > 1e6) break;
    for(int i=-2050; i < 2050; i++){
       for(int j=-2050; j < 2050; j++){
           int part = partition(i, j);
                if(part == 1) gr1->SetPoint(gr1->GetN(), i, j); 
           else if(part == 2) gr2->SetPoint(gr2->GetN(), i, j); 
           else if(part == 3) gr3->SetPoint(gr3->GetN(), i, j); 
           else if(part == 4) gr4->SetPoint(gr4->GetN(), i, j); 
           else if(part == 5) gr5->SetPoint(gr5->GetN(), i, j); 
           else if(part == 6) gr6->SetPoint(gr6->GetN(), i, j); 
           else if(part == 7) gr7->SetPoint(gr7->GetN(), i, j); 
           else if(part == 8) gr8->SetPoint(gr8->GetN(), i, j); 
           tot ++ ;
       }
    }
}
gr1->SetMarkerColor(1);
gr2->SetMarkerColor(2);
gr3->SetMarkerColor(3);
gr4->SetMarkerColor(4);
gr5->SetMarkerColor(5);
gr6->SetMarkerColor(6);
gr7->SetMarkerColor(7);
gr8->SetMarkerColor(8);

TCanvas* can = new TCanvas("can", "", 1000, 1000);
can->cd();
TH2F* h2 = new TH2F("h2","",1, -2050, 2050, 1, -2050, 2050); 
h2->Draw();
gr1->Draw("p");
gr2->Draw("p");
gr3->Draw("p");
gr4->Draw("p");
gr5->Draw("p");
gr6->Draw("p");
gr7->Draw("p");
gr8->Draw("p");
can->SaveAs("Ecal_check_cpp.png");
return 0;
}


int partition(float x, float y)
{
    float a1, b1, a2, b2 ;
    if (x>1850 && x < 2020){
        line_a_b(1850, 750  , 2020, 600 , a1, b1);
        line_a_b(1850, -1000, 2020, -850, a2, b2);
        if ( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1 ) return 1 ;
    }
    if (y < 1850 && x < 2020){
        line_a_b(750 , 1850,2020, 600, a1, b1);
        line_a_b(1000, 1850,2020, 850, a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1 ) return 2;
    }
    if (y < 2020 && y > 1850){
        line_a_b( -750, 1850 , -600, 2020 , a1, b1);
        line_a_b( 1000, 1850 , 850 , 2020 , a2, b2);
        if( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1 ) return 3;
    }
    if (y < 2020 && x > -1850){
        line_a_b(  -1850, 750  , -600 , 2020 , a1, b1);
        line_a_b(  -1850, 1000 , -850 , 2020 , a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1 ) return 4;
    }
    if (x > -2020 && x < -1850){
        line_a_b( -1850, -750 ,  -2020, -600 , a1, b1);
        line_a_b( -1850, 1000 ,  -2020,  850 , a2, b2);
        if( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1 ) return 5;
    }
    if (x > -2020 && y > -1850){
        line_a_b( -750 , -1850 , -2020, -600 , a1, b1);
        line_a_b( -1000, -1850 , -2020, -850 , a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1) return 6;
    }
    if (y > -2020 && y < -1850){
        line_a_b(  750 , -1850 ,  600, -2020 , a1, b1);
        line_a_b( -1000, -1850 , -850, -2020 , a2, b2);
        if( (a1*x + b1*y) < 1 && (a2*x+b2*y) < 1) return 7;
    }
    if (y > -2020 && x < 1850){
        line_a_b(  1850, -750  ,  600, -2020 , a1, b1);
        line_a_b(  1850, -1000 ,  850, -2020 , a2, b2);
        if( (a1*x + b1*y) > 1 && (a2*x+b2*y) < 1) return 8;
    }
    return -1;
}

void line_a_b(float x1, float y1, float x2, float y2, float& a, float& b)
{
    a = (y2-y1)/(x1*y2-x2*y1);
    b = (x1-x2)/(x1*y2-x2*y1);
}
