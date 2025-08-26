//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

Int_t plots_diff(){

  TString fileN01;
  //
  fileN01 = "./histAll.root";

  TFile *f1 = new TFile(fileN01.Data());

  TH1D *h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe = (TH1D*)f1->Get("h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe");
  TH1D *h1_fixed_chord_length_hist_reco_x_m_true_x_ctapipe = (TH1D*)f1->Get("h1_fixed_chord_length_hist_reco_x_m_true_x_ctapipe");
  TH1D *h1_reco_x_m_true = (TH1D*)f1->Get("h1_reco_x_m_true");
  TH1D *h1_reco_x_m_true_CTLearn = (TH1D*)f1->Get("h1_reco_x_m_true_CTLearn");
  //
  TH1D *h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe = (TH1D*)f1->Get("h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe");
  TH1D *h1_fixed_chord_length_hist_reco_y_m_true_y_ctapipe = (TH1D*)f1->Get("h1_fixed_chord_length_hist_reco_y_m_true_y_ctapipe");
  TH1D *h1_reco_y_m_true = (TH1D*)f1->Get("h1_reco_y_m_true");
  TH1D *h1_reco_y_m_true_CTLearn = (TH1D*)f1->Get("h1_reco_y_m_true_CTLearn");
  
  // 
  h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetLineColor(kBlack);
  h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetLineWidth(2.0);
  h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetMarkerColor(kBlack);
  // 
  h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetLineColor(kBlack);
  h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetLineWidth(2.0);
  h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetMarkerColor(kBlack);
  //
  h1_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetLineColor(kRed+2);
  h1_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetLineWidth(2.0);
  h1_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetMarkerColor(kRed+2);
  //
  h1_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetLineColor(kRed+2);
  h1_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetLineWidth(2.0);
  h1_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetMarkerColor(kRed+2);
  //
  h1_reco_x_m_true->SetLineColor(kGreen+2);
  h1_reco_y_m_true->SetLineColor(kGreen+2);
  h1_reco_x_m_true->SetLineWidth(2.0);
  h1_reco_y_m_true->SetLineWidth(2.0);
  h1_reco_x_m_true->SetMarkerColor(kGreen+2);
  h1_reco_y_m_true->SetMarkerColor(kGreen+2);
  //
  h1_reco_x_m_true_CTLearn->SetLineColor(kMagenta+2);
  h1_reco_y_m_true_CTLearn->SetLineColor(kMagenta+2);
  h1_reco_x_m_true_CTLearn->SetLineWidth(2.0);
  h1_reco_y_m_true_CTLearn->SetLineWidth(2.0);
  h1_reco_x_m_true_CTLearn->SetMarkerColor(kMagenta+2);
  h1_reco_y_m_true_CTLearn->SetMarkerColor(kMagenta+2);

  

  

  //
  //h1_03->SetLineColor(kBlue+2);
  //h1_03->SetLineWidth(3.0);
  //h1_03->SetMarkerColor(kBlue+2);
  //
  //h1_01->SetMaximum(500);
  //h1_01->Draw();
  //h1_02->Draw("sames");
  //h1_03->Draw("sames");
  //h1_01->GetXaxis()->SetTitle("r_true - r_reco");
  //h1_01->GetYaxis()->SetTitle("FADC counts");
  //  
  //
  //
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe, "Old ctapipe", "apl");
  leg->AddEntry(h1_fixed_chord_length_hist_reco_x_m_true_x_ctapipe, "Current ctapipe (closed MR.)", "apl");
  leg->AddEntry(h1_reco_x_m_true, "Dev ctapipe (draft MR.)", "apl");
  leg->AddEntry(h1_reco_y_m_true_CTLearn, "CTLearn ~6 M par.", "apl");
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,600);
  c1->Divide(2,1);
  //
  c1->cd(1);
  h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetTitle("");
  h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetMinimum(0.0);
  h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->SetMaximum(600.0);
  h1_not_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->Draw();
  h1_fixed_chord_length_hist_reco_x_m_true_x_ctapipe->Draw("samese");
  h1_reco_x_m_true->Draw("samese");
  h1_reco_x_m_true_CTLearn->Draw("samese");
  leg->Draw();  
  //
  c1->cd(2);
  h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetTitle("");
  h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetMinimum(0.0);
  h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->SetMaximum(600.0);
  h1_not_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->Draw();
  h1_fixed_chord_length_hist_reco_y_m_true_y_ctapipe->Draw("samese");
  h1_reco_y_m_true->Draw("samese");
  h1_reco_y_m_true_CTLearn->Draw("samese");
  leg->Draw();  

  return 0;
}
