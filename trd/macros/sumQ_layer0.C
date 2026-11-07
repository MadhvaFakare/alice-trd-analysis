// sumQ_layer0.C
//
// Per-track sum Q0+Q1+Q2 in TRD layer 0, electrons vs pions,
// read from the TRDPID flat tree written by Kangkan's PIDStudy task (trdpid_tree.root).
//
// Run (any ROOT 6, e.g. inside `alienv enter O2Physics/trd-pidstudy-local2`):
//   root -l -b -q 'sumQ_layer0.C("trdpid_tree.root")'            // all pT
//   root -l -b -q 'sumQ_layer0.C("trdpid_tree.root", 1.0, 2.0)'  // optional pT window in GeV/c
//
// Sanity check before using the plot: with no pT window the printout should give
//   electrons N = 145, mean ~ 76.96    pions N = 1365, mean ~ 65.78
// (the numbers from explore_tree.C). If not, stop and check the cuts.
//
// Output (in trd/plots/): sumQ_layer0.png / .pdf (linear y) and sumQ_layer0_log.png (log y)

#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

// Append the per-track layer-0 sum of all tracks in tree t that pass `cut`.
// TTree::Draw with "goff" draws nothing; it only fills an internal buffer (GetV1)
// with the value of the expression for every selected row.
static void collect(TTree* t, const char* cut, std::vector<double>& out)
{
  t->SetEstimate(t->GetEntries() + 1);   // buffer large enough for every row, otherwise it is truncated
  Long64_t n = t->Draw("fQ00+fQ10+fQ20", cut, "goff");
  double* v = t->GetV1();
  for (Long64_t i = 0; i < n; ++i) out.push_back(v[i]);
}

static void summary(const char* name, const std::vector<double>& q, double& mean, double& sd)
{
  mean = TMath::Mean(q.size(), q.data());
  sd   = TMath::StdDev(q.size(), q.data());
  printf("  %-9s N = %5zu   mean = %7.2f   std dev = %6.2f   error on mean = %5.2f\n",
         name, q.size(), mean, sd, sd / std::sqrt((double)q.size()));
}

void sumQ_layer0(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root", double ptMin = 0., double ptMax = 100.)
{
  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) { printf("Cannot open %s\n", fname); return; }

  // Layer 0 must have a real measurement in all three slices.
  // -1 is the "no usable tracklet" marker, so each slice needs > 0 (a bare "fQ00" would let -1 through).
  TString base  = Form("fQ00>0 && fQ10>0 && fQ20>0 && fPt>=%g && fPt<%g", ptMin, ptMax);
  // Kangkan's species label: 0 = electron, 1 = pion, -1 = unknown.
  // Careful: 1 is NOT pion in the ALICE o2::track::PID::ID convention (there 1 = muon, 2 = pion).
  TString cutE  = base + " && fSpecies==0";
  TString cutPi = base + " && fSpecies==1";

  // The file holds one trdpid tree per DF_* directory: loop over all of them.
  std::vector<double> qE, qPi;
  int nTrees = 0;
  TIter next(f->GetListOfKeys());
  while (TKey* key = (TKey*)next()) {
    if (!TString(key->GetName()).BeginsWith("DF_")) continue;
    TTree* t = (TTree*)f->Get(Form("%s/trdpid", key->GetName()));
    if (!t) continue;
    collect(t, cutE.Data(), qE);
    collect(t, cutPi.Data(), qPi);
    ++nTrees;
  }
  if (qE.empty() || qPi.empty()) { printf("No electrons or no pions selected, check the file/cuts\n"); return; }

  printf("\nPer-track Q0+Q1+Q2, TRD layer 0, pT in [%g, %g) GeV/c, %d trees read\n", ptMin, ptMax, nTrees);
  double mE, sE, mPi, sPi;
  summary("electrons", qE, mE, sE);
  summary("pions", qPi, mPi, sPi);

  // Same definitions as in trd_q0_comparison2.C / trd_q_slices.C, so the numbers are comparable
  double sepPower = std::fabs(mE - mPi) / (0.5 * (sE + sPi));
  double signif   = std::fabs(mE - mPi) / std::sqrt(sE * sE / qE.size() + sPi * sPi / qPi.size());
  printf("  electron excess in mean     = %5.1f %%\n", 100. * (mE - mPi) / mPi);
  printf("  separation power dMean/<SD> = %5.3f\n", sepPower);
  printf("  mean difference             = %5.1f sigma\n", signif);

  // Pion efficiency at 90% electron efficiency, unbinned:
  // threshold = value that keeps 90% of electrons above it; count pions above the same value.
  std::vector<double> sortedE(qE);
  std::sort(sortedE.begin(), sortedE.end());
  double thr = sortedE[(size_t)std::floor(0.10 * sortedE.size())];
  size_t nPiAbove = std::count_if(qPi.begin(), qPi.end(), [thr](double q) { return q >= thr; });
  double effPi = (double)nPiAbove / qPi.size();
  double effPiErr = std::sqrt(effPi * (1. - effPi) / qPi.size());
  printf("  at 90%% electron efficiency (cut Q >= %.1f): pion efficiency = %.3f +- %.3f, rejection = %.2f\n\n",
         thr, effPi, effPiErr, effPi > 0 ? 1. / effPi : 0.);

  // ---------------- histograms ----------------
  // Q0, Q1 < 127 and Q2 < 62 before correction; the corrected sum can go a bit higher, so 0-400.
  const int nb = 40; const double lo = 0., hi = 400.;
  TH1D* hE  = new TH1D("hSumQ_El", "", nb, lo, hi);
  TH1D* hPi = new TH1D("hSumQ_Pi", "", nb, lo, hi);
  hE->SetDirectory(nullptr);  hPi->SetDirectory(nullptr);   // keep them alive after the file closes
  hE->Sumw2();  hPi->Sumw2();                               // proper errors after scaling
  for (double q : qE)  hE->Fill(q);
  for (double q : qPi) hPi->Fill(q);
  if (hE->GetBinContent(nb + 1) > 0 || hPi->GetBinContent(nb + 1) > 0)
    printf("WARNING: entries above %g (overflow), widen the histogram range\n", hi);
  hE->Scale(1. / qE.size());   // unit area, so shapes can be compared despite 145 vs 1365 tracks
  hPi->Scale(1. / qPi.size());

  gStyle->SetOptStat(0);
  const int colE = kRed + 1, colPi = kBlue + 1;
  hE->SetLineColor(colE);   hE->SetMarkerColor(colE);   hE->SetMarkerStyle(20); hE->SetLineWidth(2);
  hPi->SetLineColor(colPi); hPi->SetMarkerColor(colPi); hPi->SetMarkerStyle(21); hPi->SetLineWidth(2);
  hPi->GetXaxis()->SetTitle("Q_{0}+Q_{1}+Q_{2}, TRD layer 0 (corrected, a.u.)");
  hPi->GetYaxis()->SetTitle("normalised counts");
  hPi->GetXaxis()->SetTitleSize(0.045); hPi->GetYaxis()->SetTitleSize(0.045);
  hPi->GetYaxis()->SetTitleOffset(1.3);

  for (int logy = 0; logy <= 1; ++logy) {
    TCanvas* c = new TCanvas(Form("c%d", logy), "", 900, 650);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04);
    c->SetLogy(logy);
    double ymax = 1.35 * std::max(hE->GetMaximum(), hPi->GetMaximum());
    hPi->SetMaximum(logy ? 5 * ymax : ymax);
    hPi->SetMinimum(logy ? 1e-4 : 0.);
    hPi->Draw("E1");
    hE->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.52, 0.66, 0.94, 0.86);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.037);
    leg->AddEntry(hE,  Form("e  (N = %zu, #LT Q #GT = %.1f)", qE.size(), mE), "lep");
    leg->AddEntry(hPi, Form("#pi  (N = %zu, #LT Q #GT = %.1f)", qPi.size(), mPi), "lep");
    leg->Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.036);
    tx.DrawLatex(0.17, 0.85, "PbPb 2023, work in progress");
    tx.DrawLatex(0.17, 0.80, "V0 daughters tagged by TPC n#sigma");
    if (ptMin > 0. || ptMax < 100.)
      tx.DrawLatex(0.17, 0.75, Form("%.1f #leq p_{T} < %.1f GeV/c", ptMin, ptMax));

    TString out = "/home/madhvafakare/alice/my-analysis/trd/plots/";
    out += logy ? "sumQ_layer0_log" : "sumQ_layer0";
    if (ptMin > 0. || ptMax < 100.) out += Form("_pt%.1f-%.1f", ptMin, ptMax);
    c->SaveAs(out + ".png");
    if (!logy) c->SaveAs(out + ".pdf");
  }
  f->Close();
}