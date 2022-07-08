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
  // double vertth = 0.136 / 24.8;

  // double vetoth = 0.04 / 24.8;
  
  double vertth = 0.136;
 
  double vetoth = 0.04 ;
  //std::cout<<i<<" "<<j<<std::endl;
 
  //vert cut here
  int up  = j+1;
  int down = j-1;
  if(up > 63) up -= 64;
  if(down < 0) down += 64;
  //std::cout<<"vert start"<<std::endl;
  if(e[i][up] < vertth) return false;
  if(e[i][down] < vertth) return false;
 
  //veto cut here
  //std::cout<<"veto start"<<std::endl;
  if(i != 0){
   
      if( e[i - 1][j] > vetoth) return false;
      if( e[i - 1][up] > vetoth) return false;
      if( e[i - 1][down] > vetoth) return false;
  }
  if(i != 23){
      if( e[i + 1][j] > vetoth) return false;
      if( e[i + 1][up] > vetoth) return false;
      if( e[i + 1][down] > vetoth) return false;

  }
  //std::cout<<"veto end"<<std::endl;
  return true;
}

void oHCalsimulationanalysisalltower(){
  //double calibsim  = 24.8;
  double calibsim  = 1.;
  TFile* out = new TFile("oHCalcosmiccalibresultalltower.root","RECREATE");
  
  TH1D* he = new TH1D("muon_truth_e","", 300,0,100);
  TH1D* ave = new TH1D("avg_e","", 300,0,1);
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
  TH1D* hcalibfl = new TH1D("calibfl","", 300,0,1);
  TH1D* hphi = new TH1D("phi","", 64,0,2*M_PI);
  TH1D* heta = new TH1D("eta","", 24,-1.1,1.1);
  TH1D* htheta = new TH1D("theta","", 100,0.1,M_PI / 2+0.1);
  //make tower histograms here
  TH1D* hrveto[24][64];
  for(int i = 0; i < 24; i++){
    for(int j = 0; j < 64; j++){
      std::ostringstream hvertvetoname;
      hvertvetoname<<"h_raw_vert_veto_"<<i<<"_"<<j;
      // hrveto[i][j] = new TH1D(hvertvetoname.str().c_str(),"",50,0, 1);
      hrveto[i][j] = new TH1D(hvertvetoname.str().c_str(),"",180,0, 1.065/calibsim);
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

  TH1D* energy = new TH1D("partE","", 400, 0.5, 200.5);
  energy->GetXaxis()->SetTitle("muon energy(GeV)");
  energy->GetYaxis()->SetTitle("muon/GeV/m^2/sr");
  htheta->GetXaxis()->SetTitle("theta");
  htheta->GetYaxis()->SetTitle("count");
  htower->GetXaxis()->SetTitle("ieta");
  htower->GetYaxis()->SetTitle("iphi");
  henergy2->GetXaxis()->SetTitle("towerE");
  henergy2->GetYaxis()->SetTitle("muonE");
  hcalibsim->GetXaxis()->SetTitle("tower_calib");
  hcalibsim->GetYaxis()->SetTitle("tower_sim");
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
  t->Add("/sphenix/user/bseidlitz/cosmicCalib/hcalSim/triggerAna/rootFiles/hcal_in_out_2/*.root");
  t->Add("/sphenix/user/bseidlitz/cosmicCalib/hcalSim/triggerAna/rootFiles/hcal_in_out_1/*.root");
  //t->Add("2Mtest1.root");
  //t->Add("2Mtest2.root");
  //t->Add("testtest5.root");
  //t->Add("testtest6.root");
  
  int tower_calib_n;
  int tower_sim_n;
  float part_E;
  float part_px;
  float part_py;
  float part_pz;
  float truth_vx;
  float truth_vy;
  float truth_vz;
  
  float tower_calib_E[100000];
  float tower_calib_eta[100000];
  float tower_calib_phi[100000];
  int tower_calib_iphi[100000];
  int tower_calib_ieta[100000];

  float tower_sim_E[100000];
  float tower_sim_eta[100000];
  float tower_sim_phi[100000];
  int tower_sim_iphi[100000];
  int tower_sim_ieta[100000];
  
  t->SetBranchAddress("tower_calib_n", &tower_calib_n);
  t->SetBranchAddress("tower_sim_n", &tower_sim_n);
  t->SetBranchAddress("tower_calib_E", &tower_calib_E);
  t->SetBranchAddress("tower_sim_E", &tower_sim_E);
  t->SetBranchAddress("part_E", &part_E);

  t->SetBranchAddress("part_px", &part_px);
  t->SetBranchAddress("part_py", &part_py);
  t->SetBranchAddress("part_pz", &part_pz);
  t->SetBranchAddress("truth_vx", &truth_vx);
  t->SetBranchAddress("truth_vy", &truth_vy);
  t->SetBranchAddress("truth_vz", &truth_vz);
  t->SetBranchAddress("tower_calib_eta", &tower_calib_eta);
  t->SetBranchAddress("tower_calib_phi", &tower_calib_phi);
  t->SetBranchAddress("tower_calib_ieta", &tower_calib_ieta);
  t->SetBranchAddress("tower_calib_iphi", &tower_calib_iphi);

  t->SetBranchAddress("tower_sim_eta", &tower_sim_eta);
  t->SetBranchAddress("tower_sim_phi", &tower_sim_phi);
  // t->SetBranchAddress("tower_sim_ieta", &tower_sim_ieta);
  // t->SetBranchAddress("tower_sim_iphi", &tower_sim_iphi);
  int ne = t->GetEntries();
  for(int e=0;e<ne;e++){
  //for(int e=0;e<20000;e++){
    if(e%10000==0) std::cout<<"number "<<e<<"/"<<(t->GetEntries())<<std::endl;
    t->GetEntry(e);
    double totalE[24][64] = {{0}};
    double totalEcalib[24][64] = {{0}};
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
    if(fillhis) energy->Fill(part_E, weight);
     
    for(int n = 0; n<tower_calib_n;n++){
      //totalE[tower_sim_ieta[n]][tower_sim_iphi[n]] += tower_sim_E[n];
      totalEcalib[tower_calib_ieta[n]][tower_calib_iphi[n]] += tower_calib_E[n];
  
    }
    int nactive = 0;
    bool iffill = true;
    for(int i = 0; i<24;i++){
      for(int j = 0; j<64; j++){
	if(totalEcalib[i][j]>(0.04 / calibsim)){
	  ave->Fill(totalE[i][j]);
	  if(totalEcalib[i][j]>0.04){
	    //hcalibsim->Fill(totalE[i][j], totalEcalib[i][j]);
	    if((totalE[i][j]>0.019)&&(totalE[i][j]<0.021))hcalibfl->Fill(totalEcalib[i][j]);
	  }
	    // he->Fill(totalE[i][j]);
	  htower->Fill(i,j);
	 
	  nactive++;
	  //std::cout<<"i = "<<i<<" j = "<<std::j<<std::endl;
	  if(fillmip(totalEcalib, i, j)){
	  
	    
	    hrveto[i][j]->Fill(totalEcalib[i][j]);
	    // std::cout<<"i = "<<i<<" j = "<<j<<std::endl;
	    if(i!=0 && i!=23){
	      henergy2->Fill(totalEcalib[i][j], part_E);
	      if(totalE[i][j]>0) he->Fill(part_E);
	    }
	  }
	  
	  
	}
      }
    }
 
    hna->Fill(nactive+0.5);
  }
  energy->Fit("fit","","",0.5,200);
  htheta->Fit("fittheta");
  out->Write();
  out->Close();
}
