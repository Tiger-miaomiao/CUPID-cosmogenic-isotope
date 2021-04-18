#include "TH1D.h"
#include "spectra.h"

void flux()
{

TH1D* hneutron=new TH1D("hneutron","hneutron",63,neutronX);
TH1D* hproton=new TH1D("hproton","hproton",23,protonX);
TH1D* hgamma=new TH1D("hgamma","hgamma",22,gammaX);
TH1D* hmuon=new TH1D("hmuon","hmuon",22,muonX);

 for(int i=1;i<=63;i++)
   {
     hneutron->SetBinContent(i,(neutronY[i-1]+neutronY[i])/2.0);
   }

 for(int i=1;i<=23;i++)
   {
     hproton->SetBinContent(i,(protonY[i-1]+protonY[i])/2.0);
   }

 for(int i=1;i<=22;i++)
   {
     hgamma->SetBinContent(i,(gammaY[i-1]+gammaY[i])/2.0);
   }

 for(int i=1;i<=22;i++)
   {
     hmuon->SetBinContent(i,(muonY[i-1]+muonY[i])/2.0);
   }


 //////////////////////draw plot

 TCanvas* c1=new TCanvas("c1","c1",0,0,600,400);
 c1->SetLogy();
 c1->SetLogx();

 hneutron->GetYaxis()->SetTitle("Flux (cm^{-2}s^{-1}MeV^{-1})");
 hneutron->GetXaxis()->SetTitle("Energy (MeV)");

 hneutron->Draw("");
 hproton->Draw("same");
 hgamma->Draw("same");
 hmuon->Draw("same");



//////////////////////calculate flux

double flux_neutron=0;
double flux_proton=0;
double flux_muon=0;
double flux_gamma=0;

 for(int i=39;i<=63;i++){flux_neutron=flux_neutron+(hneutron->GetBinContent(i)*hneutron->GetBinWidth(i));}
 for(int i=1;i<=23;i++){flux_proton=flux_proton+(hproton->GetBinContent(i)*hproton->GetBinWidth(i));}
 for(int i=1;i<=22;i++){flux_muon=flux_muon+(hmuon->GetBinContent(i)*hmuon->GetBinWidth(i));}
 for(int i=1;i<=22;i++){flux_gamma=flux_gamma+(hgamma->GetBinContent(i)*hgamma->GetBinWidth(i));}

 cout<<"neutron flux:"<<flux_neutron<<"/cm^{2}/s"<<endl;
 cout<<"proton flux:"<<flux_proton<<"/cm^{2}/s"<<endl;
 cout<<"muon flux:"<<flux_muon<<"/cm^{2}/s"<<endl;
 cout<<"gamma flux:"<<flux_gamma<<"/cm^{2}/s"<<endl;

}


