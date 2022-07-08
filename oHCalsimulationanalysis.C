#include <TH2D.h>
#include<TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <sstream> //std::ostringstsream
#include <fstream> //std::ifstream
#include <iostream> //std::cout, std::endl
//#include <math.h>
#include <TGraphErrors.h>
#include<TGraph.h>
#include<TSpectrum.h>
#include <vector>
#include <memory>
#include <TF1.h>
Double_t fitf(Double_t* x, Double_t * par){
 double Ec = 0.5;
 double E = 4.29;
 double eps = 854;
 double n = 3.01;
 double N = 47.12402;
 double I0 = 70.7;
 double I = par[0] * N * I0 *TMath::Power(E + x[0] , -n) /(1 + x[0] / eps);
 return I;
}

bool fillmip(double e[24][64],int i, int j){
  //double vertth = 0.136 / 24.8;
  double vertth = 320;
  //double vetoth = 0.04 / 24.8;
  double vetoth = 96;
  //std::cout<<i<<" "<<j<<std::endl;
  if(j> 31) return false;
  if(j>1&&j < 30) return false;
  //vert cut here
  //std::cout<<"vert start"<<std::endl;
  if(j == 0 && e[i][1] < vertth) return false;
  if(j == 1 && e[i][0] < vertth) return false;
  if(j == 31 && e[i][30] < vertth) return false;
  if(j == 30 && e[i][31] < vertth) return false;
  //veto cut here
  //std::cout<<"veto start"<<std::endl;
  if(i != 0){
    if(j==0 || j==1){
      if( e[i - 1][0] > vetoth) return false;
      if( e[i - 1][1] > vetoth) return false;
    }
    if(j==30 || j==31){
      if( e[i - 1][30] > vetoth) return false;
      if( e[i - 1][31] > vetoth) return false;
    }
  }
  if(i != 23){
    if(j==0 || j==1){
      if( e[i + 1][0] > vetoth) return false;
      if( e[i + 1][1] > vetoth) return false;
    }
    if(j==30 || j==31){
      if( e[i + 1][30] > vetoth) return false;
      if( e[i + 1][31] > vetoth) return false;
    }
  }
  //std::cout<<"veto end"<<std::endl;
  return true;
}

void oHCalsimulationanalysis(){
  //double calibsim  = 1.;
  double calibsim  = 24.8;
  TFile* out = new TFile("oHCalcosmicrawresult.root","RECREATE");
  
  TH1D* he = new TH1D("muon_truth_e","", 300,0,100);
  TH1D* ave = new TH1D("avg_e","", 300,0,1);
  TH1D* aveside = new TH1D("avg_e_side","", 300,0,1);
  TH1D* avenside = new TH1D("avg_e_nonside","", 300,0,1);
  
  TH1D* towerraw = new TH1D("h_raw","", 50,0,1400);
  TH1D* towervert = new TH1D("h_vert","", 50,0,1400);
  TH1D* towerveto = new TH1D("h_vert_veto","", 50,0,1500);
  towervert->SetLineColor(1);
  towerveto->SetLineColor(2);
  TH1D* te = new TH1D("active_truthe","", 100,0,10);
  TH1D* hna = new TH1D("nactive","", 10,0,10);
  TH2D* htower = new TH2D("multi","", 24,0,24,64,0,64);
  TH2D* henergy2 = new TH2D("energy2","", 50,0,1,50,0,100);
  TH2D* hcalibsim = new TH2D("calib_sim","", 300,0,0.05,300,0,1);
  TH1D* hphi = new TH1D("phi","", 64,0,2*M_PI);
  TH1D* heta = new TH1D("eta","", 24,-1.1,1.1);
  TH1D* htheta = new TH1D("theta","", 100,0.1,M_PI / 2+0.1);
  //make tower histograms here
  TH1D* hrveto[48][2];
  double binmax = 0.2e-3/3.38021e-02 * 180;
  for(int i = 0; i < 48; i++){
    for(int j = 0; j < 2; j++){
      std::ostringstream hvertvetoname;
      hvertvetoname<<"h_raw_vert_veto_"<<i<<"_"<<j;
      hrveto[i][j] = new TH1D(hvertvetoname.str().c_str(),"",320,0.5, 3200.5);
      //hrveto[i][j] = new TH1D(hvertvetoname.str().c_str(),"",200,0,0.05);
      hrveto[i][j]->GetXaxis()->SetTitle("energy");
      hrveto[i][j]->GetYaxis()->SetTitle("counts");
    }
  }

  TF1* fittheta = new TF1("fittheta","[0]*cos(x)*cos(x)",0, M_PI / 2);
  TF1* fit = new TF1("fit", fitf, 0.1,100.1, 1);

  fit->SetParameter(0, 80);
  fit->SetParNames("time");
  fit->Write();
  const int nbins = 100;
  Double_t *xbin = new Double_t[nbins];
  for(int i = 0; i< nbins; i++){
    double exp = -1 + 4./nbins *i;
    xbin[i] = TMath::Power(10,exp);
   
  }
  //TH1D* energy = new TH1D("partE","", nbins, xbin);

  TH1D* energy = new TH1D("partP","", 400, 0.5, 200.5);
  energy->GetXaxis()->SetTitle("muon momentum(GeV/c)");
  energy->GetYaxis()->SetTitle("muon/GeV/c/m^2/sr");
  htheta->GetXaxis()->SetTitle("theta");
  htheta->GetYaxis()->SetTitle("count");
  htower->GetXaxis()->SetTitle("ieta");
  htower->GetYaxis()->SetTitle("iphi");
  henergy2->GetXaxis()->SetTitle("towerE");
  henergy2->GetYaxis()->SetTitle("muonE");
  hcalibsim->GetYaxis()->SetTitle("tower_calib");
  hcalibsim->GetXaxis()->SetTitle("tower_sim");
  towerveto->GetXaxis()->SetTitle("ADC");
  towerveto->GetYaxis()->SetTitle("count/s");
  he->GetXaxis()->SetTitle("GeV");
  
  double xmin = -200;
  double xmax = 0;
  double zmin = -200;
  double zmax = 0;
  double yfix = -160;

  double A = (xmax - xmin) / 100 * (zmax - zmin) / 100; 
  double ti = 581.9;
  double weight = 2 / A /0.09546;
  std::cout<<"weight = "<<weight<<std::endl;

  TChain* t = new TChain("T");
  //t->Add("testtest1.root");
  //t->Add("testtest2.root");
  //t->Add("testtest3.root");
  //t->Add("testtest4.root");
  t->Add("/sphenix/user/shuhangli/hcalfullsim/rootFiles/hcal_in_out_1/*.root");
  //t->Add("2Mtest1.root");
  //t->Add("2Mtest2.root");
  //t->Add("/sphenix/user/bseidlitz/cosmicCalib/hcalSim/triggerAna/rootFiles/hcal_in_out_2/*.root");
  //t->Add("/sphenix/user/bseidlitz/cosmicCalib/hcalSim/triggerAna/rootFiles/hcal_in_out_1/*.root");
  //t->Add("testtest5.root");
  //t->Add("testtest6.root");
  
  int tower_raw_n;
  int tower_sim_n;
  float part_E;
  float part_px;
  float part_py;
  float part_pz;
  float truth_vx;
  float truth_vy;
  float truth_vz;
  
  float tower_raw_E[100000];
  float tower_raw_eta[100000];
  float tower_raw_phi[100000];
  int tower_raw_iphi[100000];
  int tower_raw_ieta[100000];

  float tower_sim_E[100000];
  float tower_sim_eta[100000];
  float tower_sim_phi[100000];
  int tower_sim_iphi[100000];
  int tower_sim_ieta[100000];
  
  t->SetBranchAddress("tower_raw_n", &tower_raw_n);
  t->SetBranchAddress("tower_sim_n", &tower_sim_n);
  t->SetBranchAddress("tower_raw_E", &tower_raw_E);
  t->SetBranchAddress("tower_sim_E", &tower_sim_E);
  t->SetBranchAddress("part_E", &part_E);

  t->SetBranchAddress("part_px", &part_px);
  t->SetBranchAddress("part_py", &part_py);
  t->SetBranchAddress("part_pz", &part_pz);
  t->SetBranchAddress("truth_vx", &truth_vx);
  t->SetBranchAddress("truth_vy", &truth_vy);
  t->SetBranchAddress("truth_vz", &truth_vz);
  t->SetBranchAddress("tower_raw_eta", &tower_raw_eta);
  t->SetBranchAddress("tower_raw_phi", &tower_raw_phi);
  t->SetBranchAddress("tower_raw_ieta", &tower_raw_ieta);
  t->SetBranchAddress("tower_raw_iphi", &tower_raw_iphi);

  t->SetBranchAddress("tower_sim_eta", &tower_sim_eta);
  t->SetBranchAddress("tower_sim_phi", &tower_sim_phi);
  t->SetBranchAddress("tower_sim_ieta", &tower_sim_ieta);
  t->SetBranchAddress("tower_sim_iphi", &tower_sim_iphi);
  int ne = t->GetEntries();
  for(int e=0;e<ne;e++){
  // for(int e=0;e<100;e++){
    if(e%10000==0) std::cout<<"number "<<e<<"/"<<(t->GetEntries())<<std::endl;
    t->GetEntry(e);
    double totalE[24][64] = {{0}};
    double totalEraw[24][64] = {{0}};
    //energy->Fill(part_E);
    double prho = sqrt(part_px*part_px + part_pz*part_pz);

    double theta = abs(atan(prho / part_py));
    htheta->Fill(theta, 1/sin(theta));
    bool fillhis = true;
    
    if(truth_vy < yfix) fillhis = false;
    double l = (truth_vy - yfix) / abs(part_py);
    double x = truth_vx + l * part_px;
    double z = truth_vz + l * part_pz;
    if((x> xmax)||(x<xmin)) fillhis = false;
    if((z> zmax)||(z<zmin)) fillhis = false;
    //if(fillhis) htheta->Fill(theta, 1/sin(theta));
    if(theta > M_PI / 18) fillhis = false;
    part_E = sqrt(part_px * part_px + part_py * part_py + part_pz * part_pz);
    if(fillhis) energy->Fill(part_E, weight);
     
    for(int n = 0; n<tower_raw_n;n++){
      
      totalEraw[tower_raw_ieta[n]][tower_raw_iphi[n]] += tower_raw_E[n];
  
    }
    for(int n = 0; n<tower_sim_n;n++){
      totalE[tower_sim_ieta[n]][tower_sim_iphi[n]] += tower_sim_E[n];
    }

    int nactive = 0;
    bool iffill = true;
    for(int i = 0; i<24;i++){
      for(int j = 0; j<64; j++){
	if(totalEraw[i][j]>(96)){
	  //if(totalE[i][j]>(0.04 / calibsim)){
	  ave->Fill(totalEraw[i][j]);
	  if(i==0||i==23){
	    aveside->Fill(totalEraw[i][j]);
	  }
	  if(i==1 || i==22){
	    avenside->Fill(totalEraw[i][j]);
	  }
	  if((totalEraw[i][j]>0.04)&&(totalE[i][j]>0.04/24.8)) hcalibsim->Fill(totalE[i][j], totalEraw[i][j]);
	  // he->Fill(totalE[i][j]);
	  htower->Fill(i,j);
	 
	  nactive++;
	  //std::cout<<"i = "<<i<<" j = "<<std::j<<std::endl;
	  //if(fillmip(totalEcalib, i, j)){
	  if(fillmip(totalEraw, i, j)){
	    int idx1 = 0;
	    if(j > 1) idx1 = 1;
	    int idx0 = 2 * i + abs((j-1) % 2);
	    if(j>1) idx0 = 2 * i + j % 2;
	    //std::cout<<i<<" "<<j<<" "<<idx0<<" "<<idx1<<std::endl;
	    //hrveto[idx0][idx1]->Fill(totalEcalib[i][j]);
	    hrveto[idx0][idx1]->Fill(totalEraw[i][j]);
	    if(i!=0 && i!=23){
	      henergy2->Fill(totalEraw[i][j], part_E);
	      if(totalE[i][j]>0) he->Fill(part_E);
	    }
	  }
	  
	  
	}
      }
    }
    if(totalE[11][0] < 0.004) iffill = false;
    if(iffill) towerraw->Fill(1434.55 * totalE[11][0], 1./ti);
    if(totalE[11][1] < 0.136) iffill = false;
    if(iffill) towervert->Fill(1434.55* totalE[11][0], 1./ti);
    if(totalE[10][0] > 0.04) iffill = false;
    if(totalE[12][0] > 0.04) iffill = false;
    if(totalE[12][1] > 0.04) iffill = false;
    if(totalE[10][1] > 0.04) iffill = false;
    if(iffill) towerveto->Fill(1434.55* totalE[11][0], 1./ti);
    if(nactive>=2) te->Fill(part_E);
    //if(part_E>1) he->Fill(part_E);
    hna->Fill(nactive+0.5);
  }
  energy->Fit("fit","","",0.5,200);
  htheta->Fit("fittheta");
  out->Write();
  out->Close();
}
