#include "TPaveLabel.h"
#include "TPaveText.h"
void testflux(){

  TH2D *NDFlux_293kA;
  TH2D *NDFlux_280kA;
  TH1D *FDFlux;

  TH1D *StroboFlux_bin0;
  TH1D *StroboFlux_bin1;
  TH1D *StroboFlux_bin2;
  TH1D *StroboFlux_bin3;
  TH1D *StroboFlux_bin4;
  TH1D *StroboFlux_bin5;
  TH1D *StroboFlux_bin6;
  TH1D *StroboFlux_underbin;
  TH1D *Enu_all;

  TFile *f1 = new TFile("LBNF_TDRFlux_Nov19.root");
  NDFlux_293kA = (TH2D*)gDirectory->Get("ND_293kA_nu_numu");
  NDFlux_280kA = (TH2D*)gDirectory->Get("ND_280kA_nu_numu");
 

  
  TFile *f2 = new TFile("/dune/app/users/sganguly/newdev_august19/g4lbne/example/ND_Laura_sliced_flux.root");
  StroboFlux_bin0 = (TH1D*)gDirectory->Get("BeamEspect_vmu_bin0");
  StroboFlux_bin1 = (TH1D*)gDirectory->Get("BeamEspect_vmu_bin1");
  StroboFlux_bin2 = (TH1D*)gDirectory->Get("BeamEspect_vmu_bin2");
  StroboFlux_bin3 = (TH1D*)gDirectory->Get("BeamEspect_vmu_bin3");
  StroboFlux_bin4 = (TH1D*)gDirectory->Get("BeamEspect_vmu_bin4");
  StroboFlux_bin5 = (TH1D*)gDirectory->Get("BeamEspect_vmu_bin5");
  StroboFlux_bin6 = (TH1D*)gDirectory->Get("BeamEspect_vmu_bin6");
  StroboFlux_underbin = (TH1D*)gDirectory->Get("BeamEspect_vmu_under");
  
  
  std::cout<<" ND flux energy bins = "<<NDFlux_293kA->GetXaxis()->GetNbins()<<std::endl;
  
  
  TH1D *h_293;
  h_293 = NDFlux_293kA->ProjectionX("ND_293kA", 1, 1);
  int binmax1 = h_293->GetMaximumBin();
  double val1 = h_293->GetBinContent(binmax1);
  h_293->Scale(1/val1);

  TH1D *h_280;
  h_280 = NDFlux_280kA->ProjectionX("ND_280kA", 1, 1);
  int binmax2 = h_280->GetMaximumBin();
  double val2 = h_280->GetBinContent(binmax2);
  h_280->Scale(1/val2);

  h_293->SetLineColor(kBlack);
  h_293->SetLineWidth(3);
  h_280->SetLineColor(kViolet);
  h_280->SetLineWidth(3);
  
  
  int binmax_b0 = StroboFlux_bin0->GetMaximumBin();
  double valb0 = StroboFlux_bin0->GetBinContent(binmax_b0);
  StroboFlux_bin0->Scale(1/valb0);
 
  int binmax_b1 = StroboFlux_bin1->GetMaximumBin();
  double valb1 = StroboFlux_bin1->GetBinContent(binmax_b1);
  StroboFlux_bin1->Scale(1/valb1);

  int binmax_b2 = StroboFlux_bin2->GetMaximumBin();
  double valb2 = StroboFlux_bin2->GetBinContent(binmax_b2);
  StroboFlux_bin2->Scale(1/valb2);

  int binmax_b3 = StroboFlux_bin3->GetMaximumBin();
  double valb3 = StroboFlux_bin3->GetBinContent(binmax_b3);
  StroboFlux_bin3->Scale(1/valb3);

  int binmax_b4 = StroboFlux_bin4->GetMaximumBin();
  double valb4 = StroboFlux_bin4->GetBinContent(binmax_b4);
  StroboFlux_bin4->Scale(1/valb4);

  int binmax_b5 = StroboFlux_bin5->GetMaximumBin();
  double valb5 = StroboFlux_bin5->GetBinContent(binmax_b5);
  StroboFlux_bin5->Scale(1/valb5);

  int binmax_b6 = StroboFlux_bin6->GetMaximumBin();
  double valb6 = StroboFlux_bin6->GetBinContent(binmax_b6);
  StroboFlux_bin6->Scale(1/valb6);
   
  int binmax_bunder = StroboFlux_underbin->GetMaximumBin();
  double valbunder = StroboFlux_underbin->GetBinContent(binmax_bunder);
  StroboFlux_underbin->Scale(1/valbunder);
   

  StroboFlux_bin0->SetLineColor(kBlue);
  StroboFlux_bin0->SetLineWidth(3);
  StroboFlux_bin1->SetLineColor(kBlue+2);
  StroboFlux_bin1->SetLineWidth(3);
  StroboFlux_bin2->SetLineColor(kBlue-7);
  StroboFlux_bin2->SetLineWidth(3);
  StroboFlux_bin3->SetLineColor(kAzure+1);
  StroboFlux_bin3->SetLineWidth(3);
  StroboFlux_bin4->SetLineColor(kCyan);
  StroboFlux_bin4->SetLineWidth(3);
  StroboFlux_bin5->SetLineColor(kGreen+2);
  StroboFlux_bin5->SetLineWidth(3);
  StroboFlux_bin6->SetLineColor(kGreen+6);
  StroboFlux_bin6->SetLineWidth(3);
  StroboFlux_underbin->SetLineColor(kGreen-7);
  StroboFlux_underbin->SetLineWidth(3);
  
  TCanvas* c = new TCanvas("c"," ",0.,0.,900,600);
  c->cd();
  gROOT->ForceStyle();
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();

  StroboFlux_bin0->Draw("hist p L");
  StroboFlux_bin1->Draw("SAME hist p L");
  StroboFlux_bin2->Draw("SAME hist p L");
  StroboFlux_bin3->Draw("SAME hist p L");
  StroboFlux_bin4->Draw("SAME hist p L");
  StroboFlux_bin5->Draw("SAME hist p L");
  StroboFlux_bin6->Draw("SAME hist p L");
  StroboFlux_underbin->Draw("SAME hist p L");
  
  h_293->Draw("SAME hist p L");
  h_280->Draw("SAME hist p L");      
  
}
