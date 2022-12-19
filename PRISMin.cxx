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

size_t OnAxisBin;
size_t MaxOffAxisBin;
size_t NEnergyBins;
size_t NCoeffs;

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

Eigen::VectorXd Solve(Eigen::MatrixXd const &NDFluxMatrix,
                      Eigen::MatrixXd const &RegMatrix,
                      Eigen::VectorXd const &TargetFlux,
                      std::pair<double, double> const &EnuFitLimits) {

  Eigen::MatrixXd WeightingMatrix =
      Eigen::MatrixXd::Identity(NEnergyBins, NEnergyBins);

  size_t FitRegionELow =
      NDFlux_293kA->GetXaxis()->FindFixBin(EnuFitLimits.first);
  size_t FitRegionEHigh =
      NDFlux_293kA->GetXaxis()->FindFixBin(EnuFitLimits.second);

  for (int ebin = 0; ebin < NEnergyBins; ebin++) {
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
                  Eigen::MatrixXd RegMatrix, Eigen::VectorXd const &TargetFlux,
                  std::pair<double, double> const &EnuFitLimits) {

  RegMatrix *= RegParam;

  Eigen::VectorXd OffAxisWeights =
      Solve(NDFluxMatrix, RegMatrix, TargetFlux, EnuFitLimits);

  Eigen::VectorXd SolutionFlux = NDFluxMatrix * OffAxisWeights;

  TH1D AltHCAxisWeights_h("AltHCAxisWeights_h", ";280 kA;PRISM Coefficient", 1,
                          0, 1);
  AltHCAxisWeights_h.SetBinContent(1, OffAxisWeights(0));
  TH1D OffAxisWeights_h("OffAxisWeights_h",
                        ";293 kA, Off Axis (m);PRISM Coefficient", NCoeffs - 1,
                        NDFlux_293kA->GetYaxis()->GetBinLowEdge(OnAxisBin),
                        NDFlux_293kA->GetYaxis()->GetBinUpEdge(MaxOffAxisBin));
  for (int i = 1; i < NCoeffs; ++i) {
    OffAxisWeights_h.SetBinContent(i, OffAxisWeights(i));
  }

  TH1D SolutionFlux_h("SolutionFlux_h", ";280 kA, On Axis;PRISM Coefficient",
                      NEnergyBins, NDFlux_293kA->GetXaxis()->GetBinLowEdge(1),
                      NDFlux_293kA->GetXaxis()->GetBinUpEdge(NEnergyBins));

  for (int ebin = 0; ebin < NEnergyBins; ebin++) {
    SolutionFlux_h.SetBinContent(ebin + 1, SolutionFlux(ebin));
  }

  TCanvas c2("c2", "");

  double maxc = std::max(GetMax(AltHCAxisWeights_h), GetMax(OffAxisWeights_h));
  double minc = std::min(GetMin(AltHCAxisWeights_h), GetMin(OffAxisWeights_h));

  double absmaxc = std::max(std::abs(maxc), std::abs(minc));

  TPad p1("p1", "", 0, 0, 0.3, 0.5);
  p1.AppendPad();
  p1.cd();

  p1.SetLeftMargin(0.25);
  AltHCAxisWeights_h.SetNdivisions(101);
  AltHCAxisWeights_h.GetYaxis()->SetRangeUser(-absmaxc * 1.1, absmaxc * 1.1);
  AltHCAxisWeights_h.Draw();

  c2.cd();
  TPad p2("p2", "", 0.3, 0, 1, 0.5);
  p2.AppendPad();
  p2.cd();

  p2.SetLeftMargin(0);
  OffAxisWeights_h.GetYaxis()->SetRangeUser(-absmaxc * 1.1, absmaxc * 1.1);
  OffAxisWeights_h.Draw();

  c2.cd();

  TPad p3("p3", "", 0, 0.5, 1, 1);
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

  TLine l1(EnuFitLimits.first, 0, EnuFitLimits.first, maxh);
  TLine l2(EnuFitLimits.second, 0, EnuFitLimits.second, maxh);

  l1.SetLineColor(kRed);
  l2.SetLineColor(kRed);

  l1.Draw();
  l2.Draw();

  c2.cd();

  TLegend leg(0.6, 0.6, 0.9, 0.75);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  leg.AddEntry(FDFlux, "Target Flux", "l");
  leg.AddEntry(&SolutionFlux_h, "PRISM Fit", "l");

  leg.Draw();

  c2.Print("Coeffs.pdf");
}

int main(int argc, const char *argv[]) {

  gStyle->SetOptStat(false);

  TFile fin(argv[1], "READ");

  NDFlux_293kA = fin.Get<TH2D>("ND_293kA_nu_numu");
  NDFlux_280kA = fin.Get<TH2D>("ND_280kA_nu_numu");
  FDFlux = fin.Get<TH1D>("FD_nu_numu");

  OnAxisBin = NDFlux_293kA->GetYaxis()->FindFixBin(0.0);
  MaxOffAxisBin = NDFlux_293kA->GetYaxis()->FindFixBin(28.5);

  // The fluxes we are using here are very finely binned in energy, average over
  // the bins and rescale back to the same normalization
  int EnergyRebin = 10;
  NDFlux_293kA->RebinX(EnergyRebin);
  NDFlux_293kA->Scale(1.0 / double(EnergyRebin));
  NDFlux_280kA->RebinX(EnergyRebin);
  NDFlux_280kA->Scale(1.0 / double(EnergyRebin));

  // Oscillate the FD flux
  for (int e = 0; e < FDFlux->GetXaxis()->GetNbins(); ++e) {
    double bc, be;

    bc = FDFlux->GetBinContent(e + 1);
    be = FDFlux->GetBinError(e + 1);

    bc *= Oscillate(FDFlux->GetXaxis()->GetBinCenter(e + 1));

    FDFlux->SetBinContent(e + 1, bc);
  }

  FDFlux->RebinX(EnergyRebin);
  FDFlux->Scale(1.0 / double(EnergyRebin));

  NEnergyBins = FDFlux->GetXaxis()->GetNbins();

  size_t NAdditionalFluxes = 1; // Just 280kA at this point

  // Add aditional flux components here and increment the number of additional
  // fluxes

  NCoeffs = MaxOffAxisBin - OnAxisBin + NAdditionalFluxes;

  Eigen::MatrixXd NDFluxMatrix = Eigen::MatrixXd::Zero(NEnergyBins, NCoeffs);

  Eigen::VectorXd TargetFlux = Eigen::VectorXd::Zero(NEnergyBins);

  for (int ebin = 0; ebin < NEnergyBins; ++ebin) {

    // Add the on-axis 280kA flux as the first column
    NDFluxMatrix(ebin, 0) = NDFlux_280kA->GetBinContent(ebin + 1, OnAxisBin);

    // Add additional flux components to the flux matrix here
    // NDFluxMatrix(ebin, 1) = SomeFlux->GetBinContent(ebin + 1);

    TargetFlux(ebin) = FDFlux->GetBinContent(ebin + 1);

    for (int oabin = OnAxisBin; oabin < MaxOffAxisBin; ++oabin) {
      NDFluxMatrix(ebin, NAdditionalFluxes + oabin - OnAxisBin) =
          NDFlux_293kA->GetBinContent(ebin + 1, oabin);
    }
  }

  Eigen::MatrixXd RegMatrix = Eigen::MatrixXd::Zero(NCoeffs, NCoeffs);

  RegMatrix(0, 0) = 1;
  // Add additional flux component regularisation here
  //  RegMatrix(1, 1) = 1;

  for (int coeff = NAdditionalFluxes; coeff < NCoeffs - 1; ++coeff) {
    RegMatrix(coeff, coeff) = 1;
    RegMatrix(coeff, coeff + 1) = -1;
  }
  RegMatrix(NCoeffs - 1, NCoeffs - 1) = 1;

  TCanvas c1("c1", "");
  c1.Print("Coeffs.pdf[");
  SolveAndDraw(NDFluxMatrix, 1E-7, RegMatrix, TargetFlux, {0.8, 3.5});
  SolveAndDraw(NDFluxMatrix, 1E-8, RegMatrix, TargetFlux, {0.8, 3.5});
  SolveAndDraw(NDFluxMatrix, 1E-9, RegMatrix, TargetFlux, {0.8, 3.5});
  SolveAndDraw(NDFluxMatrix, 1E-10, RegMatrix, TargetFlux, {0.8, 3.5});
  c1.Print("Coeffs.pdf]");
}