#include "Eigen/Dense"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

#include "BargerPropagator.h"

#include <iostream>
#include <memory>
#include "TMath.h"

enum NuTypes {
  kNuebarType = -1,
  kNumubarType = -2,
  kNutaubarType = -3,
  kNueType = 1,
  kNumuType = 2,
  kNutauType = 3,
};

double Oscillate(double ENu_GeV, int FromType = kNumuType,
                 int ToType = kNumuType) {
  double baseline_km = 1300;

  // Sin^2(Theta_12)
  double s2th12 = 0.297;
  // Sin^2(Theta_13)
  double s2th13 = 0.0214;
  // Sin^2(Theta_23)
  double s2th23 = 0.526;
  // Dm^2_21
  double dm2_21 = 7.37E-5;
  //|Dm^2_Atm|
  double dm2_atm = 2.463E-3;

  double dcp = 0;

  double params[] = {s2th12, s2th13, s2th23, dm2_21, dm2_atm, dcp};

  static const double deg2rad = asin(1) / 90.0;
  const static double REarth_km = 6371.393;
  double DipAngle = asin(baseline_km / (2.0 * REarth_km)) / deg2rad;
  double LengthParam = cos((90.0 + DipAngle) * deg2rad);
  BargerPropagator bp;
  bp.DefinePath(LengthParam, 0);

  bp.SetMNS(s2th12, s2th13, s2th23, dm2_21, dm2_atm, dcp, ENu_GeV, true,
            FromType);
  bp.propagate(ToType);

  return bp.GetProb(FromType, ToType);
}

TH2D *NDFlux_293kA;
TH2D *NDFlux_280kA;
TH1D *FDFlux;

std::array<TH1D *, 15> NDFluxOnAxisStroboTimeBin;

struct BinMapping {
  size_t NCoeffs;
  size_t NNDAltHCBins;
  size_t NNDStroboTimeBins;
  size_t NNDOffAxisBins;
  size_t NEnergyBins;
};

double GetMax(TH1 const &h) {
  double max = -std::numeric_limits<double>::max();
  for (int i = 0; i < h.GetXaxis()->GetNbins(); ++i) {
    max = std::max(max, h.GetBinContent(i + 1));
  }
  return max;
}

double GetMin(TH1 const &h) {
  double min = std::numeric_limits<double>::max();
  for (int i = 0; i < h.GetXaxis()->GetNbins(); ++i) {
    min = std::min(min, h.GetBinContent(i + 1));
  }
  return min;
}

double GetMax(std::vector<TH1 *> const &hs) {
  double max = -std::numeric_limits<double>::max();
  for (auto const &h : hs) {
    max = std::max(max, GetMax(*h));
  }
  return max;
}
double GetMin(std::vector<TH1 *> const &hs) {
  double min = std::numeric_limits<double>::max();
  for (auto const &h : hs) {
    min = std::min(min, GetMin(*h));
  }
  return min;
}

Eigen::VectorXd Solve(Eigen::MatrixXd const &NDFluxMatrix,
                      Eigen::MatrixXd const &RegMatrix,
                      Eigen::VectorXd const &TargetFlux, BinMapping const &bm,
                      std::pair<double, double> const &EnuFitLimits) {

  Eigen::MatrixXd WeightingMatrix =
    Eigen::MatrixXd::Identity(bm.NEnergyBins, bm.NEnergyBins);

  size_t FitRegionELow =
   NDFlux_293kA->GetXaxis()->FindFixBin(EnuFitLimits.first);
  size_t FitRegionEHigh =
    NDFlux_293kA->GetXaxis()->FindFixBin(EnuFitLimits.second);

  for (int ebin = 0; ebin < bm.NEnergyBins; ebin++) {
    if (ebin <= FitRegionELow) { // low energy bin(s) weight
      WeightingMatrix(ebin, ebin) *= 0.8;
    }
    if (ebin >= FitRegionEHigh) { // high energy bin(s) weight
      WeightingMatrix(ebin, ebin) *= 0.0;
    }
  }

  return ((NDFluxMatrix.transpose() * WeightingMatrix * NDFluxMatrix) +
          RegMatrix.transpose() * RegMatrix)
    .inverse() *
    NDFluxMatrix.transpose() * WeightingMatrix * TargetFlux;
}

void SolveAndDraw(Eigen::MatrixXd const &NDFluxMatrix, double RegParam,
                  Eigen::VectorXd const &TargetFlux, BinMapping const &bm,
                  std::pair<double, double> const &EnuFitLimits) {

  Eigen::MatrixXd RegMatrix = Eigen::MatrixXd::Zero(bm.NCoeffs, bm.NCoeffs);

  for (int coeff = 0; coeff < (bm.NNDAltHCBins + bm.NNDStroboTimeBins);
       ++coeff) {
    RegMatrix(coeff, coeff) = 1;
  }

  for (int coeff = (bm.NNDAltHCBins + bm.NNDStroboTimeBins);
       coeff < bm.NCoeffs - 1; ++coeff) {
    RegMatrix(coeff, coeff) = 1;
    RegMatrix(coeff, coeff + 1) = -1;
  }
  RegMatrix(bm.NCoeffs - 1, bm.NCoeffs - 1) = 1;

  RegMatrix *= RegParam;

  Eigen::VectorXd OffAxisWeights =
    Solve(NDFluxMatrix, RegMatrix, TargetFlux, bm, EnuFitLimits);

  Eigen::VectorXd SolutionFlux = NDFluxMatrix * OffAxisWeights;

  TH1D CoeffWeights_h("CoeffWeights_h", ";C_{i};Weight", bm.NCoeffs, 0,
                      bm.NCoeffs);
  for (int i = 0; i < bm.NCoeffs; ++i) {
    CoeffWeights_h.SetBinContent(i, OffAxisWeights(i));
  }

  TH1D SolutionFlux_h("SolutionFlux_h", ";280 kA, On Axis;PRISM Coefficient",
                      bm.NEnergyBins,
                      NDFlux_293kA->GetXaxis()->GetBinLowEdge(1),
                      NDFlux_293kA->GetXaxis()->GetBinUpEdge(bm.NEnergyBins));

  for (int ebin = 0; ebin < bm.NEnergyBins; ebin++) {
    SolutionFlux_h.SetBinContent(ebin + 1, SolutionFlux(ebin));
  }

  TCanvas c2("c2", "");

  double maxc = GetMax(CoeffWeights_h);
  double minc = GetMin(CoeffWeights_h);

  double absmaxc = std::max(std::abs(maxc), std::abs(minc));

  TPad p1("p1", "", 0, 0, 1, 0.5);
  p1.AppendPad();
  p1.cd();

  CoeffWeights_h.SetNdivisions(101);
  CoeffWeights_h.SetLineWidth(2);
  CoeffWeights_h.GetYaxis()->SetRangeUser(-absmaxc * 1.1, absmaxc * 1.1);
  CoeffWeights_h.Draw();

  if (bm.NNDAltHCBins) {
    TLine *l = new TLine(bm.NNDAltHCBins, -absmaxc * 1.1, bm.NNDAltHCBins,
                         absmaxc * 1.1);
    l->SetLineColor(kRed);
    l->Draw();
  }

  if (bm.NNDStroboTimeBins) {
    TLine *l = new TLine(bm.NNDAltHCBins + bm.NNDStroboTimeBins, -absmaxc * 1.1,
                         bm.NNDAltHCBins + bm.NNDStroboTimeBins, absmaxc * 1.1);
    l->SetLineColor(kRed);
    l->Draw();
  }

  c2.cd();

  TPad p3("p3", "", 0, 0.5, 1, 1.0);
  p3.AppendPad();
  p3.cd();

  double maxh = std::max(GetMax(*FDFlux), GetMax(SolutionFlux_h));

  FDFlux->GetYaxis()->SetRangeUser(0, maxh * 1.1);

  FDFlux->SetLineWidth(2);
  FDFlux->Draw("HIST");
  SolutionFlux_h.SetLineColor(kGreen);
  SolutionFlux_h.Draw("HIST SAME");
  

  TLatex tlx;
  tlx.SetTextAlign(32);

  std::stringstream ss("");
  ss << "Regularization Strength: " << std::scientific << RegParam;

  tlx.DrawLatexNDC(0.89, 0.85, ss.str().c_str());
  tlx.DrawLatexNDC(0.89, 0.8,
                   ("NAltHCBins: " + std::to_string(bm.NNDAltHCBins)).c_str());
  tlx.DrawLatexNDC(
      0.89, 0.75,
      ("NStroboBins: " + std::to_string(bm.NNDStroboTimeBins)).c_str());
  tlx.DrawLatexNDC(
      0.89, 0.7, ("NPRISMBins: " + std::to_string(bm.NNDOffAxisBins)).c_str());

  TLine l1(EnuFitLimits.first, 0, EnuFitLimits.first, maxh);
  TLine l2(EnuFitLimits.second, 0, EnuFitLimits.second, maxh);

  l1.SetLineColor(kRed);
  l2.SetLineColor(kRed);

  l1.Draw();
  l2.Draw();

  c2.cd();

  TLegend leg(0.6, 0.7, 0.9, 0.85);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  leg.AddEntry(FDFlux, "Target Flux", "l");
  leg.AddEntry(&SolutionFlux_h, "PRISM Fit", "l");

  leg.Draw();


  /*TPad p6("p6", "", 0, 0.3, 1, 0.6);
  p6.AppendPad();
  p6.cd();
  TH1F *ratio = (TH1F*)SolutionFlux_h.Clone("ratio");
  ratio->Divide(FDFlux);
  ratio->SetLineWidth(2);
  ratio->SetLineColor(kBlack);
  ratio->GetYaxis()->SetTitle("PRISM Fit/Target Flux");
  ratio->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
  ratio->Draw("HIST SAME");
  */
  c2.cd();

  c2.Print("Coeffs.pdf");
}

void DrawFluxes(int OnAxisBin, int FirstFullOffAxisBin, int MaxOffAxisBin,
                int NNDStroboTimeBins) {

  TCanvas c3("c3", "");

  TH1D *reference_flux = NDFlux_293kA->ProjectionX(
      "df_reference_flux", FirstFullOffAxisBin, FirstFullOffAxisBin);

  std::cout << "Ref bin: " << FirstFullOffAxisBin << std::endl;

  // top left -- ND Off Axis + AltHC
  TPad p1("c3_p1", "", 0, 0.5, 0.5, 1);
  p1.AppendPad();
  p1.cd();

  TH1D *althc = NDFlux_280kA->ProjectionX("df_althc", 1, 1);
  std::vector<TH1 *> off_axisfluxes;
  std::vector<TH1 *> off_axisfluxes_dcs;
  for (int i = OnAxisBin; i < MaxOffAxisBin; i += 10) {
    off_axisfluxes.push_back(NDFlux_293kA->ProjectionX(
        ("df_prism_slice" + std::to_string(i)).c_str(), i, i));
  }
  for (int i = 0; i < off_axisfluxes.size(); i++) {
    if ((i == OnAxisBin)) {
      double yaxup = 1.1 * GetMax(off_axisfluxes);
      double yaxlow = 0.9 * GetMin(off_axisfluxes);
      off_axisfluxes[i]->GetYaxis()->SetRangeUser(yaxlow, yaxup);
      off_axisfluxes[i]->GetYaxis()->SetTitle("Flux (A.U.)");
      off_axisfluxes_dcs.push_back(
          static_cast<TH1 *>(off_axisfluxes[i]->DrawClone("CHIST PLC")));
    } else {
      off_axisfluxes_dcs.push_back(
          static_cast<TH1 *>(off_axisfluxes[i]->DrawClone("CHIST PLC SAME")));
    }
  }
  reference_flux->SetLineColor(kBlack);
  reference_flux->DrawClone("CHIST SAME");
  althc->SetLineColor(kRed);
  althc->DrawClone("CHIST SAME");

  TLegend *leg1 = new TLegend(0.4, 0.6, 0.95, 0.9);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);

  leg1->AddEntry(reference_flux, "Reference flux (2.25m off axis)", "l");
  leg1->AddEntry(althc, "280 kA", "l");
  leg1->AddEntry(off_axisfluxes_dcs[2], "293 kA, 10m off axis", "l");

  leg1->Draw();

  c3.cd();
  // top right -- Strobo
  TPad p2("c3_p2", "", 0.5, 0.5, 1, 1);
  p2.AppendPad();
  p2.cd();

  std::vector<TH1 *> strobo_fluxes;
  std::vector<TH1 *> strobo_fluxes_dcs;
  for (int i = 0; i < NNDStroboTimeBins; i++) {
    strobo_fluxes.push_back(
        static_cast<TH1 *>(NDFluxOnAxisStroboTimeBin[i]->Clone(
            ("df_strobo_slice" + std::to_string(i)).c_str())));
  }

  for (int i = 0; i < NNDStroboTimeBins; i++) {
    if (!i) {
      double yaxup = 1.1 * GetMax(strobo_fluxes);
      double yaxlow = 0.9 * GetMin(strobo_fluxes);
      strobo_fluxes[i]->GetYaxis()->SetTitle("Flux (A.U.)");
      strobo_fluxes[i]->GetYaxis()->SetRangeUser(yaxlow, yaxup);
      strobo_fluxes_dcs.push_back(
				  static_cast<TH1 *>(strobo_fluxes[i]->DrawClone("CHIST PLC")));
    } else {
      strobo_fluxes_dcs.push_back(
				  static_cast<TH1 *>(strobo_fluxes[i]->DrawClone("CHIST PLC SAME")));
    }
  }

  TLegend *leg2 = new TLegend(0.4, 0.6, 0.95, 0.9);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);

  leg2->AddEntry(strobo_fluxes_dcs[0], "ND On axis Time Underflow", "l");
  leg2->AddEntry(strobo_fluxes_dcs[1], "ND On axis Time bin 2", "l");
  leg2->AddEntry(strobo_fluxes_dcs[3], "ND On axis Time bin 4", "l");

  leg2->Draw();

  c3.cd();
  // bottom left
  TPad p3("c3_p3", "", 0, 0, 0.5, 0.5);
  p3.AppendPad();
  p3.cd();
  for (int i = 0; i < off_axisfluxes.size(); i++) {
    off_axisfluxes[i]->Divide(reference_flux);
  }
  for (int i = 0; i < off_axisfluxes.size(); i++) {
    if ((i == OnAxisBin)) {
      double yaxup = 1.1 * GetMax(off_axisfluxes);
      double yaxlow = 0.9 * GetMin(off_axisfluxes);
      off_axisfluxes[i]->GetYaxis()->SetTitle("Flux/Reference");
      off_axisfluxes[i]->GetYaxis()->SetRangeUser(yaxlow, yaxup);
      off_axisfluxes[i]->DrawClone("CHIST PLC");
    } else {
      off_axisfluxes[i]->DrawClone("CHIST PLC SAME");
    }
  }
  althc->Divide(reference_flux);
  althc->DrawClone("CHIST SAME");

  c3.cd();
  // bottom right
  TPad p4("c3_p4", "", 0.5, 0, 1, 0.5);
  p4.AppendPad();
  p4.cd();
  for (int i = 0; i < NNDStroboTimeBins; i++) {
    strobo_fluxes[i]->Divide(reference_flux);
  }
  for (int i = 0; i < NNDStroboTimeBins; i++) {
    if (!i) {
      double yaxup = 1.1 * GetMax(strobo_fluxes);
      double yaxlow = 0.9 * GetMin(strobo_fluxes);
      strobo_fluxes[i]->GetYaxis()->SetTitle("Flux/Reference");
      strobo_fluxes[i]->GetYaxis()->SetRangeUser(yaxlow, yaxup);
      strobo_fluxes[i]->DrawClone("CHIST PLC");
    } else {
      strobo_fluxes[i]->DrawClone("CHIST PLC SAME");
    }
  }

  c3.Print("InputFluxes.pdf");
}

int main(int argc, const char *argv[]) {

  gStyle->SetOptStat(false);

  TFile PRISMFluxesFile(argv[1], "READ");
  TFile StroboNDFluxesFile(argv[2], "READ");

  NDFlux_293kA = PRISMFluxesFile.Get<TH2D>("ND_293kA_nu_numu");
  NDFlux_280kA = PRISMFluxesFile.Get<TH2D>("ND_280kA_nu_numu");
  FDFlux = PRISMFluxesFile.Get<TH1D>("FD_nu_numu");
  //FDFlux = PRISMFluxesFile.Get<TH1D>("BeamEspect_vmu_bin0");

  size_t NAltHCBins = 1;
  size_t OnAxisBin = NDFlux_293kA->GetYaxis()->FindFixBin(0.0);
  size_t FirstFullOffAxisBin = NDFlux_293kA->GetYaxis()->FindFixBin(2.1);
  size_t MaxOffAxisBin = NDFlux_293kA->GetYaxis()->FindFixBin(28.5);
  size_t NPRISMBins = MaxOffAxisBin - OnAxisBin;
  size_t NPRISMOffAxisBins = MaxOffAxisBin - FirstFullOffAxisBin;

  // 14 additional histograms from strobocopic approach 
  size_t NNDStroboTimeBins = 15;
  NDFluxOnAxisStroboTimeBin[0] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_under");
  NDFluxOnAxisStroboTimeBin[1] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin0");
  NDFluxOnAxisStroboTimeBin[2] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin1");
  NDFluxOnAxisStroboTimeBin[3] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin2");
  NDFluxOnAxisStroboTimeBin[4] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin3");
  NDFluxOnAxisStroboTimeBin[5] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin4");
  NDFluxOnAxisStroboTimeBin[6] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin5");
  NDFluxOnAxisStroboTimeBin[7] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin6");
  NDFluxOnAxisStroboTimeBin[8] =
  StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin7");
  NDFluxOnAxisStroboTimeBin[9] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin8");
  NDFluxOnAxisStroboTimeBin[10] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin9");
  NDFluxOnAxisStroboTimeBin[11] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin10");
  NDFluxOnAxisStroboTimeBin[12] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin11");
  NDFluxOnAxisStroboTimeBin[13] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin12");
  NDFluxOnAxisStroboTimeBin[14] =
    StroboNDFluxesFile.Get<TH1D>("BeamEspect_vmu_bin13");
  
  // The fluxes we are using here are very finely binned in energy, average over
  // the bins and rescale back to the same normalization
  int EnergyRebin = 1;
  NDFlux_293kA->RebinX(EnergyRebin);
  NDFlux_293kA->Scale(1.0 / double(EnergyRebin));
  NDFlux_280kA->RebinX(EnergyRebin);
  NDFlux_280kA->Scale(1.0 / double(EnergyRebin));

  // Rebin the alt HC run to average over the whole on axis position
  int AltHCOffAxisRebin = 4;
  NDFlux_280kA->RebinY(AltHCOffAxisRebin);
  NDFlux_280kA->Scale(1.0 / double(AltHCOffAxisRebin));

  int EnergyRebinStrobo = EnergyRebin;
  for (int i = 0; i < NNDStroboTimeBins; ++i) {
    NDFluxOnAxisStroboTimeBin[i]->SetTitle("");
    NDFluxOnAxisStroboTimeBin[i]->RebinX(EnergyRebinStrobo);
    NDFluxOnAxisStroboTimeBin[i]->Scale(1.0 / double(EnergyRebinStrobo));
  }

  // Oscillate the FD flux
  for (int e = 0; e < FDFlux->GetXaxis()->GetNbins(); ++e) {
    double bc, be;

    bc = FDFlux->GetBinContent(e + 1);
    be = FDFlux->GetBinError(e + 1);

    bc *= Oscillate(FDFlux->GetXaxis()->GetBinCenter(e + 1), kNumubarType,
                    kNumubarType);

    FDFlux->SetBinContent(e + 1, bc);
  }

  FDFlux->RebinX(EnergyRebin);
  FDFlux->Scale(1.0 / double(EnergyRebin));

  FDFlux->Reset();

  double GaussMean = 1.5; //GeV

  //Overwrite the oscillated flux with a pure gaussian
  for(int i = 0; i < FDFlux->GetXaxis()->GetNbins(); ++i){
    double Enu = FDFlux->GetXaxis()->GetBinCenter(i+1);
    FDFlux->SetBinContent(i+1, TMath::Gaus(Enu, GaussMean, GaussMean*0.15, true));
  }



  // End loading and rebinning of fluxes
  DrawFluxes(OnAxisBin, FirstFullOffAxisBin, MaxOffAxisBin, NNDStroboTimeBins);

  // peak scale hte FD and ND strobo fluxes
  FDFlux->Scale(1.0 / GetMax(*FDFlux));
  for (int i = 0; i < NNDStroboTimeBins; ++i) {
    NDFluxOnAxisStroboTimeBin[i]->Scale(1.0 /
                                        GetMax(*NDFluxOnAxisStroboTimeBin[i]));
  }

  size_t NEnergyBins = FDFlux->GetXaxis()->GetNbins();

  // make matrices and bin mappings for a few different configurations
  Eigen::MatrixXd NDFluxMatrix_PRISMAltHC =
    Eigen::MatrixXd::Zero(NEnergyBins, NPRISMBins + NAltHCBins);
  BinMapping NDBinMapping_PRISMAltHC{NPRISMBins + NAltHCBins, NAltHCBins, 0,
      NPRISMBins, NEnergyBins};

  Eigen::MatrixXd NDFluxMatrix_PRISMNoAltHC =
    Eigen::MatrixXd::Zero(NEnergyBins, NPRISMBins);
  BinMapping NDBinMapping_PRISMNoAltHC{NPRISMBins, 0, 0, NPRISMBins,
      NEnergyBins};

  Eigen::MatrixXd NDFluxMatrix_PRISMStroboAltHC = Eigen::MatrixXd::Zero(
									NEnergyBins, NPRISMOffAxisBins + NAltHCBins + NNDStroboTimeBins);
  BinMapping NDBinMapping_PRISMStroboAltHC{
    NPRISMOffAxisBins + NAltHCBins + NNDStroboTimeBins, NAltHCBins,
      NNDStroboTimeBins, NPRISMOffAxisBins, NEnergyBins};

  Eigen::MatrixXd NDFluxMatrix_PRISMStroboNoAltHC =
    Eigen::MatrixXd::Zero(NEnergyBins, NPRISMOffAxisBins + NNDStroboTimeBins);
  BinMapping NDBinMapping_PRISMStroboNoAltHC{
    NPRISMOffAxisBins + NNDStroboTimeBins, 0, NNDStroboTimeBins,
      NPRISMOffAxisBins, NEnergyBins};

  Eigen::MatrixXd NDFluxMatrix_NoPRISMStroboAltHC =
    Eigen::MatrixXd::Zero(NEnergyBins, NNDStroboTimeBins + NAltHCBins);
  BinMapping NDBinMapping_NoPRISMStroboAltHC{NNDStroboTimeBins + NAltHCBins, NAltHCBins, NNDStroboTimeBins, 0,
      NEnergyBins};

  Eigen::VectorXd TargetFlux = Eigen::VectorXd::Zero(NEnergyBins);

  // peak normalize the 2D histograms
  double AltHCNormFact = -std::numeric_limits<double>::max();
  for (int ebin = 0; ebin < NEnergyBins; ++ebin) {
    AltHCNormFact =
      std::max(AltHCNormFact, NDFlux_280kA->GetBinContent(ebin + 1, 1));
  }
  std::map<int, double> OffAxisNormFact;
  for (int oabin = OnAxisBin + 1; oabin < MaxOffAxisBin; ++oabin) {
    OffAxisNormFact[oabin] = -std::numeric_limits<double>::max();
    for (int ebin = 0; ebin < NEnergyBins; ++ebin) {
      OffAxisNormFact[oabin] = std::max(
					OffAxisNormFact[oabin], NDFlux_293kA->GetBinContent(ebin + 1, oabin));
    }
  }

  // fill the ND matrices
  for (int ebin = 0; ebin < NEnergyBins; ++ebin) {

    // each matrix tracks which bins it has filled as it goes along
    int PRISMAltHC_fluxbin = 0;
    int PRISMNoAltHC_fluxbin = 0;
    int PRISMStroboAltHC_fluxbin = 0;
    int PRISMStroboNoAltHC_fluxbin = 0;
    int NoPRISMStroboAltHC_fluxbin = 0;

    // Add the on-axis 280kA flux as the first column
    double bc = NDFlux_280kA->GetBinContent(ebin + 1, 1) / AltHCNormFact;
    // for matrices that include AltHC, add them first
    NDFluxMatrix_PRISMAltHC(ebin, PRISMAltHC_fluxbin++) = bc;
    NDFluxMatrix_PRISMStroboAltHC(ebin, PRISMStroboAltHC_fluxbin++) = bc;
    NDFluxMatrix_NoPRISMStroboAltHC(ebin, NoPRISMStroboAltHC_fluxbin++) = bc;

    for (size_t StroboNDTimeBin = 0; StroboNDTimeBin < NNDStroboTimeBins;
         ++StroboNDTimeBin) {
      double bc =
	NDFluxOnAxisStroboTimeBin[StroboNDTimeBin]->GetBinContent(ebin + 1);
      // for matrices that include Strobo, add them next
      NDFluxMatrix_PRISMStroboAltHC(ebin, PRISMStroboAltHC_fluxbin++) = bc;
      NDFluxMatrix_PRISMStroboNoAltHC(ebin, PRISMStroboNoAltHC_fluxbin++) = bc;
      NDFluxMatrix_NoPRISMStroboAltHC(ebin, NoPRISMStroboAltHC_fluxbin++) = bc;
    }

    TargetFlux(ebin) = FDFlux->GetBinContent(ebin + 1);

    for (int oabin = OnAxisBin + 1; oabin < MaxOffAxisBin; ++oabin) {
      double bc =
	NDFlux_293kA->GetBinContent(ebin + 1, oabin) / OffAxisNormFact[oabin];

      // Add the PRISM off axis slices

      NDFluxMatrix_PRISMAltHC(ebin, PRISMAltHC_fluxbin++) = bc;
      NDFluxMatrix_PRISMNoAltHC(ebin, PRISMNoAltHC_fluxbin++) = bc;

      // For strobo we want to ignore the full on axis prediction
      if (oabin >= FirstFullOffAxisBin) {
        NDFluxMatrix_PRISMStroboAltHC(ebin, PRISMStroboAltHC_fluxbin++) = bc;
        NDFluxMatrix_PRISMStroboNoAltHC(ebin, PRISMStroboNoAltHC_fluxbin++) =
	  bc;
	//NDFluxMatrix_NoPRISMStroboAltHC(ebin, NoPRISMStroboAltHC_fluxbin++) =
	//bc;
      }
    }

    if (!ebin) {
      std::cout << "PRISMAltHC_fluxbin: " << PRISMAltHC_fluxbin << std::endl;
      std::cout << "PRISMNoAltHC_fluxbin: " << PRISMNoAltHC_fluxbin
                << std::endl;
      std::cout << "PRISMStroboAltHC_fluxbin: " << PRISMStroboAltHC_fluxbin
                << std::endl;
      std::cout << "PRISMStroboNoAltHC_fluxbin: " << PRISMStroboNoAltHC_fluxbin
                << std::endl;
      std::cout << "NoPRISMStroboAltHC_fluxbin: " << NoPRISMStroboAltHC_fluxbin
                << std::endl;
    }
  }

  TCanvas c1("c1", "");
  c1.Print("Coeffs.pdf[");
  // SolveAndDraw(NDFluxMatrix_PRISMAltHC, 1E-2, TargetFlux,
  //              NDBinMapping_PRISMAltHC, {0.8, 5});
  // SolveAndDraw(NDFluxMatrix_PRISMAltHC, 2.5E-2, TargetFlux,
  //              NDBinMapping_PRISMAltHC, {0.8, 5});
  SolveAndDraw(NDFluxMatrix_PRISMAltHC, 5.0E-3, TargetFlux,
                NDBinMapping_PRISMAltHC, {0.6, 8.0});
  // SolveAndDraw(NDFluxMatrix_PRISMAltHC, 7.5E-2, TargetFlux,
  //              NDBinMapping_PRISMAltHC, {0.8, 5});

  // SolveAndDraw(NDFluxMatrix_PRISMNoAltHC, 1E-2, TargetFlux,
  //              NDBinMapping_PRISMNoAltHC, {0.8, 5});
  SolveAndDraw(NDFluxMatrix_PRISMNoAltHC, 5.0E-3, TargetFlux,
                NDBinMapping_PRISMNoAltHC, {0.6,8.0});
     //SolveAndDraw(NDFluxMatrix_PRISMNoAltHC, 5.0E-3, TargetFlux,
     //        NDBinMapping_PRISMNoAltHC, {0.8, 5.0});
  // SolveAndDraw(NDFluxMatrix_PRISMNoAltHC, 7.5E-2, TargetFlux,
  //              NDBinMapping_PRISMNoAltHC, {0.8, 5});

   SolveAndDraw(NDFluxMatrix_PRISMStroboNoAltHC, 5.0E-3, TargetFlux,
                NDBinMapping_PRISMStroboNoAltHC, {0.6, 8.0});
  //SolveAndDraw(NDFluxMatrix_PRISMStroboNoAltHC, 2.5E-2, TargetFlux,
  //              NDBinMapping_PRISMStroboNoAltHC, {0.8, 5});
  //SolveAndDraw(NDFluxMatrix_PRISMStroboNoAltHC, 5.0E-3, TargetFlux,
  //             NDBinMapping_PRISMStroboNoAltHC, {0.8, 5.0});
     //SolveAndDraw(NDFluxMatrix_PRISMStroboNoAltHC, 0.7E-2, TargetFlux,
     //         NDBinMapping_PRISMStroboNoAltHC, {0.8, 5});

   SolveAndDraw(NDFluxMatrix_PRISMStroboAltHC, 5.0E-3, TargetFlux,
                NDBinMapping_PRISMStroboAltHC, {0.6, 8.0});
  // SolveAndDraw(NDFluxMatrix_PRISMStroboAltHC, 2.5E-2, TargetFlux,
  //              NDBinMapping_PRISMStroboAltHC, {0.8, 5});
  //SolveAndDraw(NDFluxMatrix_PRISMStroboAltHC, 5.0E-3, TargetFlux,
  //             NDBinMapping_PRISMStroboAltHC, {0.8, 5.0});
     //SolveAndDraw(NDFluxMatrix_PRISMStroboAltHC, 1.0E-3, TargetFlux,
     //         NDBinMapping_PRISMStroboAltHC, {0.8, 5});

  //SolveAndDraw(NDFluxMatrix_NoPRISMStroboAltHC, 1.0E-2, TargetFlux,                                                      
  //              NDBinMapping_NoPRISMStroboAltHC, {1.0, 5.0});   
  // SolveAndDraw(NDFluxMatrix_NoPRISMStroboAltHC, 2.5E-2, TargetFlux,                                                      
  //              NDBinMapping_NoPRISMStroboNoAltHC, {0.8, 5});     
  SolveAndDraw(NDFluxMatrix_NoPRISMStroboAltHC,5E-3, TargetFlux,
               NDBinMapping_NoPRISMStroboAltHC, {0.6, 8.0});
  //SolveAndDraw(NDFluxMatrix_NoPRISMStroboAltHC, 7.5E-2, TargetFlux,                                                      
  //           NDBinMapping_NoPRISMStroboAltHC, {0.8, 5});     
   
  c1.Print("Coeffs.pdf]");
}
