#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include<TMultiGraph.h>
#include "TROOT.h"
#include<TLegend.h>
#include<TFile.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "/sphenix/u/shuhang98/AtlasStyle.C"
#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"

void comparesimulation(){
  SetAtlasStyle();

  bool dorate = 0;
  TFile* fc = new TFile("oHCalcosmicrawresult.root");
  TFile* fs = new TFile("oHCalcosmicsimresult.root");
  //TFile* fr =  new TFile("11626result.root");
  TFile* fr =  new TFile("11183result.root");
  TFile* out = new TFile("simulationcompareshape.root","RECREATE");
  TCanvas* can[48];

  for(int i = 0; i<48;i++){
    std::ostringstream cname;
    cname<<"tower_"<<i;
    std::ostringstream hname;
    hname<<"h_raw_vert_veto_"<<i<<"_0";
    std::ostringstream hrname;
    hrname<<"h_raw_vert_veto_"<<i;
    can[i] = new TCanvas(cname.str().c_str(),"",1);
    can[i]->cd();
    TH1D* hc = (TH1D*) (fc->Get(hname.str().c_str()));
    TH1D* hs = (TH1D*) (fs->Get(hname.str().c_str()));
    TH1D* hr = (TH1D*) (fr->Get(hrname.str().c_str()));
    hc->Rebin(2);
    //hs->Rebin(2);
    
    
    double scalexc = (hr->GetMean())/(hc->GetMean());
    double scalexs = (hr->GetMean())/(hs->GetMean())+100;

    double cxmax  = hc->GetXaxis()->GetXmax();
    double sxmax  = hs->GetXaxis()->GetXmax();
    int binc = hc->GetXaxis()->GetNbins();
    int bins = hs->GetXaxis()->GetNbins();
    int binr = hr->GetXaxis()->GetNbins();

    double wc = hc->GetXaxis()->GetBinWidth(1);
    double ws = hs->GetXaxis()->GetBinWidth(1);
    double wr = hr->GetXaxis()->GetBinWidth(1);
    
    hc->GetXaxis()->Set(binc, 0, scalexc * cxmax);
    hs->GetXaxis()->Set(bins, 0, scalexs * sxmax);
    
   

    hc->SetLineColor(2);
    hc->SetMarkerColor(2);
    hc->SetMarkerSize(0.7);
    hs->SetLineColor(4);
    hs->SetMarkerColor(4);
    hs->SetMarkerSize(0.7);
    hs->SetFillColor(4);
    hs->SetFillStyle(3001);
    hr->SetLineColor(1);
    hr->SetMarkerStyle(20);
    hr->SetMarkerColor(1);
    hr->SetMarkerSize(0.7);
    
    if(dorate){
      hc->Scale(1./(8437 * wc* scalexc));
      //hc->Scale(1./(583 * wc* scalexc));
      hs->Scale(1./(8437 * ws* scalexs));
      
      hr->Scale(1./(3600 * wr));
      
      hr->GetYaxis()->SetTitle("count/s/ADC");
      hs->GetYaxis()->SetTitle("count/s/ADC");
      hc->GetYaxis()->SetTitle("count/s/ADC");
    }
    else{
      hc->Scale(1./(hc->GetEntries())/(wc * scalexc));
      hs->Scale(1./(hs->GetEntries())/(ws * scalexs));
      
      hr->Scale(1./(hr->GetEntries())/wr);
      
      hr->GetYaxis()->SetTitle("#frac{dN}{NdE}");
      hs->GetYaxis()->SetTitle("#frac{dN}{NdE}");
      hc->GetYaxis()->SetTitle("#frac{dN}{NdE}");
    }
    
    hr->GetXaxis()->SetRangeUser(0,2100);
    hr->GetYaxis()->SetTitleOffset(1.8);
    //gPad->SetLeftMargin(10);
    //hs->GetYaxis()->SetRangeUser(0,0.005);
    hr->Draw("EX0");
    hs->Draw("hist same");
    //hc->Draw("EX0 same");
    //gPad->SetLogy();
    gPad->Update();
    double ypos = can[i]->GetFrame()->GetY2();
    // ypos = pow(10,ypos);
    //std::cout<<ypos<<std::endl;
    TGaxis * A = new TGaxis(0,ypos,2100,ypos,0,2100/scalexs,510,"-");
    A->SetName("axis1");
    A->SetTitle("GeV");
    A->Draw(); 
    auto legend3 = new TLegend(0.1,0.7,0.3,0.9);
    //legend3->AddEntry(hc,"simulation: tower_raw","pl");
    legend3->AddEntry(hs,"GEANT4 truth energy deposition","lf");
    legend3->AddEntry(hr,"teststand result","pl");
    legend3->Draw();
    out->cd();
    can[i]->Write();
    
  }
  out->Write();
  out->Close();
}
