//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <time.h>

//root
#include "TROOT.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1D.h"

using namespace std;

Double_t func(Double_t phi, Double_t *par);
Double_t dcor_phi0( Double_t rho, Double_t R, Double_t phi0, Double_t phi);
Double_t dcor( Double_t rho, Double_t R, Double_t phi);
void fcn(int &npar, double *gin, double &f, double *par, int iflag);

void get_initial_rho_and_phi0( Double_t R_mirror, Double_t &rho, Double_t &phi0, Double_t &phi0_min);

void fit_phi_dist_with_Minuit(Double_t amplitude_in, Double_t R_mirror_in, Double_t R_camera_in, Double_t rho_in, Double_t phi0_in, Double_t pedestal_in,
			      Double_t &amplitude_out, Double_t &R_mirror_out, Double_t &R_camera_out, Double_t &rho_out, Double_t &phi0_out, Double_t &pedestal_out,
			      Double_t &amplitude_err, Double_t &R_mirror_err, Double_t &R_camera_err, Double_t &rho_err, Double_t &phi0_err, Double_t &pedestal_err);


void read_data(TString fname, TGraphErrors *gr);
void read_file_list(TString fname, vector<TString> &file_list_v);

Double_t *fit_and_plot(TString csf_file, TCanvas *c1);

TCanvas *test_func(TString nameTitle, Double_t amplitude, Double_t R_mirror, Double_t R_camera, Double_t rho, Double_t phi0, Double_t pedestal);

TGraphErrors *_gr = new TGraphErrors();
Int_t _verbose = 0;
Int_t _event_id = 0;


int main(int argc, char *argv[]){
  //
  clock_t start, finish;
  start = clock();
  if(argc == 4 && atoi(argv[1]) == 0){
    TString inListFile = argv[2];
    TString histOut = argv[3];
    //
    cout<<"inListFile "<<inListFile<<endl
	<<"histOut    "<<histOut<<endl;
    //
    vector<TString> file_list_v;
    read_file_list(inListFile, file_list_v);
    //
    Double_t amplitude;
    Double_t R_mirror;
    Double_t R_camera;
    Double_t rho;
    Double_t phi0;
    Double_t pedestal;
    //
    TH1D *h1_amplitude = new TH1D("h1_amplitude","h1_amplitude",100, 0.0, 30.0);
    TH1D *h1_R_mirror = new TH1D("h1_R_mirror","h1_R_mirror",100, 0.0, 30.0);
    TH1D *h1_R_camera = new TH1D("h1_R_camera","h1_R_camera",100, 0.0, 10.0);
    TH1D *h1_rho = new TH1D("h1_rho","h1_rho",100, -30.0, 30.0);
    TH1D *h1_phi0 = new TH1D("h1_phi0","h1_phi0",50, -190.0, -190);
    TH1D *h1_pedestal = new TH1D("h1_pedestal","h1_pedestal",50, -10.0, 10.0);
    //
    for(unsigned int i = 0; i < file_list_v.size(); i++){
      gStyle->SetPalette(1);
      gStyle->SetFrameBorderMode(0);
      gROOT->ForceStyle();
      gStyle->SetStatColor(kWhite);
      //gStyle->SetOptStat(kFALSE);
      //c1->Divide(3,2);
      TString canvaNT = "c1";
      canvaNT += file_list_v.at(i);
      TCanvas *c1 = new TCanvas(canvaNT.Data(),canvaNT.Data(),10,10,1000,1000);
      c1->SetRightMargin(0.03);
      c1->SetLeftMargin(0.12);
      c1->SetTopMargin(0.07);
      c1->SetBottomMargin(0.1);
      c1->SetGridx();
      c1->SetGridy();
      Double_t *par_out = fit_and_plot(file_list_v.at(i),c1);
      delete c1;
      amplitude = par_out[0];
      R_mirror = par_out[1];
      R_camera = par_out[2];
      rho = par_out[3];
      phi0 = par_out[4];
      pedestal = par_out[5];
      //
      //      
      if(_verbose>1){
	cout<<"amplitude "<<amplitude<<endl
	    <<"R_mirror  "<<R_mirror<<endl
	    <<"R_camera  "<<R_camera<<endl
	    <<"rho       "<<rho<<endl
	    <<"phi0      "<<phi0<<endl
	    <<"pedestal  "<<pedestal<<endl;
      }
      //
      //
      h1_amplitude->Fill(amplitude);
      h1_R_mirror->Fill(R_mirror);
      h1_R_camera->Fill(R_camera);
      h1_rho->Fill(rho);
      h1_phi0->Fill(phi0);
      h1_pedestal->Fill(pedestal);
    }
    //
    //
    //
    TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
    rootFile->cd();
    if (rootFile->IsZombie()){
      cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
      assert(0);
    }
    else
      cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
    //
    //
    //
    //
    //test_func(TString nameTitle, Double_t amplitude, Double_t R_mirror, Double_t R_camera, Double_t rho, Double_t phi0, Double_t pedestal)
    test_func("c1_test_func_0m", 8, 11, 1.7, 0, 0.0, 0.0)->Write();
    test_func("c1_test_func_1m", 8, 11, 1.7, 1, 0.0, 0.0)->Write();
    test_func("c1_test_func_2m", 8, 11, 1.7, 2, 0.0, 0.0)->Write();
    test_func("c1_test_func_3m", 8, 11, 1.7, 3, 0.0, 0.0)->Write();
    test_func("c1_test_func_4m", 8, 11, 1.7, 4, 0.0, 0.0)->Write();
    //
    h1_amplitude->Write();
    h1_R_mirror->Write();
    h1_R_camera->Write();
    h1_rho->Write();
    h1_phi0->Write();
    h1_pedestal->Write();
    //
    rootFile->Close();
  }
  else{
    cout<<"  runID [1] = 0        "<<endl
	<<"        [2] - inListFile "<<endl
    	<<"        [2] - histOut "<<endl;
  }  //
  finish = clock();
  cout<<"-------------------------"<<endl
      <<"Working time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<endl
      <<"-------------------------"<<endl;  
  return 0;  
}

TCanvas *test_func(TString nameTitle, Double_t amplitude, Double_t R_mirror, Double_t R_camera, Double_t rho, Double_t phi0, Double_t pedestal){
  //
  Int_t nn = 1000;
  //Double_t phi_min = -180.0;
  //Double_t phi_max = 180.0;
  Double_t phi_min = -TMath::Pi();
  Double_t phi_max = TMath::Pi();
  Double_t phi;
  //
  Double_t *par = new Double_t[6];
  par[0] = amplitude;
  par[1] = R_mirror;
  par[2] = R_camera;
  par[3] = rho;
  par[4] = phi0;
  par[5] = pedestal;
  //
  Double_t *par_mirr = new Double_t[6];
  par_mirr[0] = amplitude;
  par_mirr[1] = R_mirror;
  par_mirr[2] = 0.0;
  par_mirr[3] = rho;
  par_mirr[4] = phi0;
  par_mirr[5] = pedestal;
  //
  Double_t *par_cam = new Double_t[6];
  par_cam[0] = -amplitude;
  par_cam[1] = 0.0;
  par_cam[2] = R_camera;
  par_cam[3] = rho;
  par_cam[4] = phi0;
  par_cam[5] = pedestal;
  //
  //
  TGraph *gr = new TGraph();
  TGraph *gr_mirr = new TGraph();
  TGraph *gr_cam = new TGraph();
  for(Int_t i = 0; i<nn;i++){
    phi = phi_min + (phi_max - phi_min)/(nn-1)*i;
    gr->SetPoint(gr->GetN(),
		 phi,
		 func( phi, par));
    gr_mirr->SetPoint(gr_mirr->GetN(),
		      phi,
		      func( phi, par_mirr));
    gr_cam->SetPoint(gr_cam->GetN(),
		     phi,
		     func( phi, par_cam));
  }
  //
  //
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(kBlack);
  gr_mirr->SetMarkerStyle(7);
  gr_mirr->SetMarkerColor(kRed+2);
  gr_cam->SetMarkerStyle(7);
  gr_cam->SetMarkerColor(kBlue+2);

  
  //
  //
  TCanvas *c1 = new TCanvas(nameTitle.Data(),nameTitle.Data(),10,10,1000,1000);    
  c1->SetRightMargin(0.03);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.07);
  c1->SetBottomMargin(0.1);
  c1->SetGridx();
  c1->SetGridy();
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(nameTitle.Data());
  mg->Add(gr);
  mg->Add(gr_mirr);
  mg->Add(gr_cam);
  //
  mg->GetXaxis()->SetTitle("phi, rad");
  mg->GetYaxis()->SetTitle("Integrated pixel intensity, p.e.");
  //
  mg->SetMinimum(-10.0);
  mg->SetMaximum(220.0);
  mg->Draw("APL");
  //
  return c1;
}

Double_t *fit_and_plot(TString csf_file, TCanvas *c1){
  TGraph *gr_cord = new TGraph();
  TGraph *gr_camera = new TGraph();
  gr_cord->SetNameTitle("gr_cord","gr_cord");

  Int_t nn = 1000;
  //Double_t phi_min = -180.0;
  //Double_t phi_max = 180.0;
  Double_t phi_min = -TMath::Pi();
  Double_t phi_max = TMath::Pi();
  Double_t phi;

  Double_t amplitude = 10;
  Double_t R_mirror = 11.0;
  Double_t R_camera = 1.7;
  Double_t R_camera_zero = 0;
  Double_t rho = 8.0;
  Double_t phi0 = 100.0/180.0*TMath::Pi();
  Double_t phi0_min;
  Double_t pedestal = 0.0;
  
  Double_t par[6];
  par[0] = amplitude;
  par[1] = R_mirror;
  par[2] = R_camera;
  par[3] = rho;
  par[4] = phi0;
  par[5] = pedestal;
  
  //read_data("./phi_dist/hist_phi_csvName1.csv", _gr);
  read_data(csf_file.Data(), _gr);
  //
  //
  get_initial_rho_and_phi0( R_mirror, rho, phi0, phi0_min);
  //
  //
  
  Double_t amplitude_zero_out, R_mirror_zero_out, R_camera_zero_out, rho_zero_out, phi0_zero_out, pedestal_zero_out; 
  Double_t amplitude_zero_err, R_mirror_zero_err, R_camera_zero_err, rho_zero_err, phi0_zero_err, pedestal_zero_err;
  
  fit_phi_dist_with_Minuit(amplitude, R_mirror, R_camera_zero, rho, phi0, pedestal,
			   amplitude_zero_out, R_mirror_zero_out, R_camera_zero_out, rho_zero_out, phi0_zero_out, pedestal_zero_out,
			   amplitude_zero_err, R_mirror_zero_err, R_camera_zero_err, rho_zero_err, phi0_zero_err, pedestal_zero_err);
  
  
  Double_t par_reco_zero[6];
  par_reco_zero[0] = amplitude_zero_out;
  par_reco_zero[1] = R_mirror_zero_out;
  par_reco_zero[2] = R_camera_zero_out;
  par_reco_zero[3] = rho_zero_out;
  par_reco_zero[4] = phi0_zero_out;
  par_reco_zero[5] = pedestal_zero_out;
  
  Double_t amplitude_out, R_mirror_out, R_camera_out, rho_out, phi0_out, pedestal_out; 
  Double_t amplitude_err, R_mirror_err, R_camera_err, rho_err, phi0_err, pedestal_err; 
  
  fit_phi_dist_with_Minuit(amplitude_zero_out, R_mirror, R_camera, rho_zero_out, phi0_zero_out, pedestal_zero_out,
  			   amplitude_out, R_mirror_out, R_camera_out, rho_out, phi0_out, pedestal_out,
  			   amplitude_err, R_mirror_err, R_camera_err, rho_err, phi0_err, pedestal_err);


  Double_t *par_reco = new Double_t[6];
  par_reco[0] = amplitude_out;
  par_reco[1] = R_mirror_out;
  par_reco[2] = R_camera_out;
  par_reco[3] = rho_out;
  par_reco[4] = phi0_out;
  par_reco[5] = pedestal_out;

  Double_t *par_reco_camera = new Double_t[6];
  par_reco_camera[0] = -amplitude_out;
  par_reco_camera[1] = 0.0;
  par_reco_camera[2] = R_camera_out;
  par_reco_camera[3] = rho_out;
  par_reco_camera[4] = phi0_out;
  par_reco_camera[5] = pedestal_out;

  
  
  TGraph *gr_reco = new TGraph();
  TGraph *gr_reco_zero = new TGraph();
  
  for(Int_t i = 0; i<nn;i++){
    phi = phi_min + (phi_max - phi_min)/(nn-1)*i;
    gr_cord->SetPoint(gr_cord->GetN(),
		      phi,
		      func( phi, par));
    gr_reco->SetPoint(gr_reco->GetN(),
		      phi,
		      func( phi, par_reco));
    gr_reco_zero->SetPoint(gr_reco_zero->GetN(),
			   phi,
			   func( phi, par_reco_zero));
    gr_camera->SetPoint(gr_camera->GetN(),
			phi,
			func( phi, par_reco_camera));
  }

  
  
  //gr_reco_ring->SetMarkerStyle(7);
  gr_cord->SetMarkerStyle(7);
  gr_cord->SetMarkerColor(kMagenta+2);
  _gr->SetMarkerStyle(20);
  _gr->SetMarkerColor(kBlack);
  _gr->SetLineColor(kBlack);
  _gr->SetLineWidth(2);
  gr_reco->SetMarkerStyle(7);
  gr_reco->SetMarkerColor(kRed+2);
  gr_reco_zero->SetMarkerStyle(7);
  gr_reco_zero->SetMarkerColor(kBlue+2);
  //
  //
  //
  gr_camera->SetMarkerStyle(7);
  gr_camera->SetMarkerColor(kGreen+2);
  //
  //
  //
  //
  //
  TString mg_title = csf_file;
  mg_title += " ev_id : ";
  mg_title += _event_id;
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(mg_title.Data());
  mg->Add(gr_camera);
  //mg->Add(gr_cord);
  mg->Add(gr_reco_zero);
  mg->Add(gr_reco);
  mg->Add(_gr);


  mg->GetXaxis()->SetTitle("phi, rad");
  mg->GetYaxis()->SetTitle("Integrated pixel intensity, p.e.");


  mg->SetMinimum(-10.0);
  mg->SetMaximum(220.0);
  mg->Draw("APL");

  //TLine ln_phi0(phi0,0.0,phi0,210);
  //ln_phi0.SetLineWidth(2.0);
  //ln_phi0.SetLineColor(kMagenta+2);
  //ln_phi0.Draw("same");

  //TLine ln_phi0_min(phi0_min,0.0,phi0_min,210);
  //ln_phi0_min.SetLineWidth(2.0);
  //ln_phi0_min.SetLineColor(kGreen+2);
  //ln_phi0_min.Draw("same");

  TString pdf_out = csf_file;
  pdf_out += ".pdf";
  c1->SaveAs(pdf_out.Data());
  
  return par_reco;
}

Double_t func(Double_t phi, Double_t *par){
  Double_t amplitude = par[0];
  Double_t R_mirror = par[1];
  Double_t R_camera = par[2];
  Double_t rho = par[3];
  Double_t phi0 = par[4];
  Double_t pedestal = par[5];
  //Double_t pedestal = 0.0;
  return amplitude*(dcor_phi0( rho, R_mirror, phi0, phi) - dcor_phi0( rho, R_camera, phi0, phi)) + pedestal;
}

Double_t dcor_phi0( Double_t rho, Double_t R, Double_t phi0, Double_t phi){
  Double_t d = dcor( rho, R, (phi - phi0));
  Double_t d_phi = dcor( rho, R, (phi - phi0 + 2*TMath::Pi()));
  Double_t d_phi_m = dcor( rho, R, (phi - phi0 - 2*TMath::Pi()));
  
  if(d > 0)
   return d;
  if(d_phi > 0)
    return d_phi;
  if(d_phi_m > 0)
    return d_phi_m;
  return 0.0;
}


Double_t dcor( Double_t rho, Double_t R, Double_t phi){
  if(R <= 0.0)
    return 0;
  Double_t rho_R = TMath::Abs(rho/R);
  Double_t sqrt_sin_2 = 0.0;
  Double_t omega = 1.0;
  Double_t d = 0;
  Double_t norm = 1.0;
  Double_t phi_critical;
  if(rho_R > 1.0){
    sqrt_sin_2 = rho_R*rho_R*TMath::Power(TMath::Sin(omega*(phi)),2.0);
    phi_critical = TMath::ASin(1.0/rho_R);
    if(sqrt_sin_2<=1.0){
      if(TMath::Abs(phi) < phi_critical)
	d = 2.0*norm*R*TMath::Sqrt(1.0 - sqrt_sin_2); 
      else
	d = 0.0;
    }
    else{
      d = 0.0;
    }
  }
  else{
    sqrt_sin_2 = rho_R*rho_R*TMath::Power(TMath::Sin(omega*(phi)),2.0);
    if(sqrt_sin_2<=1.0){
      d = norm*R*TMath::Sqrt(1.0 - sqrt_sin_2);
      d += norm*R*rho_R*TMath::Cos(omega*(phi));
    }
    else{
      d = 0.0;
    }
  }
  return d;
}


//event_id x y
//4600 -3.0368728984701336 52.384903
//...
void read_data(TString fname, TGraphErrors *gr){
  string mot;
  ifstream fFile(fname);
  Double_t event_id;
  Double_t x, y;
  Int_t point_counter = 0;
  if(fFile.is_open()){
    fFile>>mot>>mot>>mot;
    while(fFile>>event_id>>x>>y){
      gr->SetPoint(point_counter,x,y);
      gr->SetPointError(point_counter,6.0*TMath::Pi()/180.0,TMath::Sqrt(y));
      point_counter++;
    }
    fFile.close();
  }
  _event_id = event_id;
}

void fit_phi_dist_with_Minuit(Double_t amplitude_in, Double_t R_mirror_in, Double_t R_camera_in, Double_t rho_in, Double_t phi0_in, Double_t pedestal_in,
			      Double_t &amplitude_out, Double_t &R_mirror_out, Double_t &R_camera_out, Double_t &rho_out, Double_t &phi0_out, Double_t &pedestal_out,
			      Double_t &amplitude_err, Double_t &R_mirror_err, Double_t &R_camera_err, Double_t &rho_err, Double_t &phi0_err, Double_t &pedestal_err){
  //
  Int_t npar = 6;
  TMinuit *gMinuit = new TMinuit(npar);
  gMinuit->SetPrintLevel(-1.0);
  gMinuit->SetFCN(fcn); 
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // 
  // Set starting values and step sizes for parameters
  gMinuit->mnparm(0, "A", amplitude_in, 0.001, 0,0,ierflg);
  gMinuit->mnparm(1, "R_m", R_mirror_in, 0.0, 0,0,ierflg);
  gMinuit->mnparm(2, "R_c", R_camera_in, 0.0, 0,0,ierflg);
  gMinuit->mnparm(3, "rho", rho_in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(4, "phi0", phi0_in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(5, "pedestal", pedestal_in, 0.0, 0,0,ierflg);
  //  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  // Print results
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  //
  gMinuit->GetParameter(0, amplitude_out, amplitude_err);
  gMinuit->GetParameter(1, R_mirror_out, R_mirror_err);
  gMinuit->GetParameter(2, R_camera_out, R_camera_err);
  gMinuit->GetParameter(3, rho_out, rho_err);
  gMinuit->GetParameter(4, phi0_out, phi0_err);
  rho_out = TMath::Abs(rho_out);
  gMinuit->GetParameter(5, pedestal_out, pedestal_err);
  //
  //
  //
  //
  if(_verbose>2){
    cout<<"amplitude_out "<<amplitude_out<<endl
	<<"R_mirror_out  "<<R_mirror_out<<endl
	<<"R_camera_out  "<<R_camera_out<<endl 
	<<"rho_out       "<<rho_out<<endl
	<<"phi0_out      "<<phi0_out<<endl
	<<"pedestal_out  "<<pedestal_out<<endl;
  }
  //
  //
  //
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag){
   double chisq = 0;
   double x, y;
   double x_err, y_err;
   double tot_err;
   double delta;
   for (int i = 0; i<_gr->GetN(); i++){
     _gr->GetPoint( i, x, y);
     x_err = _gr->GetErrorX(i);
     y_err = _gr->GetErrorY(i);
     tot_err = 1.0*TMath::Sqrt(x_err*x_err + y_err*y_err);
     if(y>0){
       delta = (y - func(x, par))/tot_err;
       delta *= delta;
       chisq += delta;
     }
   }
   f = chisq;
}

void get_initial_rho_and_phi0( Double_t R_mirror, Double_t &rho, Double_t &phi0, Double_t &phi0_min){
  rho = 6.0;
  phi0 = 0.0;
  double x, y;
  double y_max_y, y_max_x;
  double y_min_y, y_min_x; 
  for (int i = 0; i<_gr->GetN(); i++){
    _gr->GetPoint( i, x, y);
    if(i == 0){
      y_max_y = y;
      y_min_y = y; 
      //
      y_max_x = x;
      y_min_x = x;
    }
    else{
      if(y_max_y<y){
	y_max_y = y;
	y_max_x = x;
      }
      if(y_min_y>y){
	y_min_y = y;
	y_min_x = x;
      }
    }
  }
  rho = (y_max_x - y_min_x)/2.0/8.0;
  phi0 = y_max_x; 
  phi0_min = y_min_x;
}

void read_file_list(TString fname, vector<TString> &file_list_v){
  string mot;
  ifstream fFile(fname);
  if(fFile.is_open()){
    while(fFile>>mot)
      file_list_v.push_back(mot);
    fFile.close();
  }
}
