// trd_q_separation.C
//
// Electron vs pion separation for ALL TRD charge slices (Q0, Q1, Q2) and all 6 layers,
// read from the flat tree written by the PIDStudy task (trdpid_tree.root).
//
// Run:   root -l -b -q 'trd_q_separation.C("trdpid_tree.root")'
// With a momentum window (fPt in GeV/c):
//        root -l -b -q 'trd_q_separation.C("trdpid_tree.root", 1.0, 2.0)'
// rej@90 = -1 in the tables means no pion survived the cut (rejection is a lower limit)
//
// Tree layout (per DF_* directory, tree "trdpid"):
//   fQ00..fQ05  slice 0 (Q0), layers 0..5     (-1 means: no usable tracklet in that layer)
//   fQ10..fQ15  slice 1 (Q1), layers 0..5
//   fQ20..fQ25  slice 2 (Q2), layers 0..5
//   fPt         track pT
//   fSpecies    0 = electron, 1 = pion, -1 = unknown

#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// Two numbers that describe how well two distributions separate.
// Assumes electrons sit at HIGHER charge than pions (transition radiation).
//   sep   = (mean_e - mean_pi) / sqrt(0.5 * (rms_e^2 + rms_pi^2))
//   rej90 = 1 / (pion efficiency) at the cut that keeps 90% of the electrons
// ---------------------------------------------------------------------------
void metrics(TH1D* e, TH1D* p, double& sep, double& rej90)
{
  sep = 0;
  rej90 = 0;
  if (e->GetEntries() < 2 || p->GetEntries() < 2) return;

  double rmsAvg = std::sqrt(0.5 * (e->GetRMS() * e->GetRMS() + p->GetRMS() * p->GetRMS()));
  if (rmsAvg > 0) sep = (e->GetMean() - p->GetMean()) / rmsAvg;

  int nb = e->GetNbinsX();
  double totE = e->Integral(1, nb);
  double totP = p->Integral(1, nb);
  if (totE <= 0 || totP <= 0) return;

  // walk the cut down from the high-charge end until 90% of electrons are kept
  double cum = 0;
  int cutBin = nb;
  for (int b = nb; b >= 1; --b) {
    cum += e->GetBinContent(b);
    if (cum / totE >= 0.90) { cutBin = b; break; }
  }
  double pionEff = p->Integral(cutBin, nb) / totP;
  if (pionEff > 0) rej90 = 1.0 / pionEff;
  else rej90 = -1; // no pion survives the cut: rejection is only a lower limit, flagged as -1
}

void trd_q_separation(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root",
                      double ptMin = 0.0, double ptMax = 1e9)
{
  gStyle->SetOptStat(0);
  gSystem->mkdir("plots", true);

  const int nS = 3; // slices Q0, Q1, Q2
  const int nL = 6; // layers

  // ----- histograms -----
  // per (slice, layer): every tracklet with q > 0
  TH1D* hE[nS][nL];
  TH1D* hP[nS][nL];
  // per track: average over the valid layers of that track (needs >= 4 valid layers)
  TH1D* hAvgE[nS + 1];
  TH1D* hAvgP[nS + 1];
  const char* sliceName[nS + 1] = {"Q0", "Q1", "Q2", "Q0+Q1+Q2"};

  for (int s = 0; s < nS; ++s) {
    for (int l = 0; l < nL; ++l) {
      hE[s][l] = new TH1D(Form("hE_%d_%d", s, l), Form("%s layer %d;charge;fraction", sliceName[s], l), 64, 0, 128);
      hP[s][l] = new TH1D(Form("hP_%d_%d", s, l), "", 64, 0, 128);
      hE[s][l]->Sumw2();
      hP[s][l]->Sumw2();
    }
  }
  for (int s = 0; s <= nS; ++s) {
    // max possible: Q0,Q1 < 127 (7-bit), Q2 < 62 (6-bit) -> sum < 316
    double hi = (s < nS) ? 128 : 320;
    hAvgE[s] = new TH1D(Form("hAvgE_%d", s), Form("%s, mean over valid layers;charge;fraction", sliceName[s]), 64, 0, hi);
    hAvgP[s] = new TH1D(Form("hAvgP_%d", s), "", 64, 0, hi);
    hAvgE[s]->Sumw2();
    hAvgP[s]->Sumw2();
  }

  // ----- open file, loop over all DF_* directories -----
  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot open %s\n", fname);
    return;
  }

  long nTracks = 0, nE = 0, nP = 0, nEUsed = 0, nPUsed = 0;

  TIter next(f->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    if (!TString(key->GetName()).BeginsWith("DF_")) continue;
    TDirectory* dir = (TDirectory*)key->ReadObj();
    TTree* t = (TTree*)dir->Get("trdpid");
    if (!t) continue;

    // leaves are read type-agnostically with GetValue(), so float/int/uchar all work
    TLeaf* q[nS][nL];
    for (int s = 0; s < nS; ++s)
      for (int l = 0; l < nL; ++l)
        q[s][l] = t->GetLeaf(Form("fQ%d%d", s, l));
    TLeaf* lPt = t->GetLeaf("fPt");
    TLeaf* lSp = t->GetLeaf("fSpecies");
    if (!lPt || !lSp || !q[0][0]) {
      printf("Missing branches in %s, skipping\n", key->GetName());
      continue;
    }

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      ++nTracks;

      int sp = (int)lSp->GetValue();
      double pt = lPt->GetValue();
      if (sp != 0 && sp != 1) continue;        // only tagged electrons and pions
      if (pt < ptMin || pt >= ptMax) continue; // optional momentum window
      if (sp == 0) ++nE; else ++nP;

      // per-layer histograms, and sums for the per-track average
      double sum[nS + 1] = {0, 0, 0, 0};
      int cnt[nS + 1] = {0, 0, 0, 0};
      int nValid = 0; // layers where all three slices are > 0

      for (int l = 0; l < nL; ++l) {
        double v[nS];
        bool all = true;
        for (int s = 0; s < nS; ++s) {
          v[s] = q[s][l]->GetValue();
          if (v[s] > 0) {
            (sp == 0 ? hE[s][l] : hP[s][l])->Fill(v[s]);
            sum[s] += v[s];
            ++cnt[s];
          } else {
            all = false;
          }
        }
        if (all) {
          sum[nS] += v[0] + v[1] + v[2];
          ++cnt[nS];
          ++nValid;
        }
      }

      // per-track average, only for tracks with >= 4 valid layers (same idea as the task's quality gate)
      for (int s = 0; s <= nS; ++s) {
        if (cnt[s] >= 4) {
          (sp == 0 ? hAvgE[s] : hAvgP[s])->Fill(sum[s] / cnt[s]);
          if (s == nS) { if (sp == 0) ++nEUsed; else ++nPUsed; }
        }
      }
    }
  }

  printf("\nTracks read: %ld | tagged electrons %ld, pions %ld (pT window %.2f to %.2f)\n",
         nTracks, nE, nP, ptMin, ptMax);
  printf("Tracks with >= 4 layers where Q0,Q1,Q2 are all > 0: electrons %ld, pions %ld\n\n", nEUsed, nPUsed);

  // ----- per-layer table -----
  printf("PER LAYER (tracklets with q > 0)\n");
  printf("%-6s %-6s %6s %6s %8s %8s %8s %8s\n", "slice", "layer", "N_e", "N_pi", "mean_e", "mean_pi", "sep", "rej@90");
  for (int s = 0; s < nS; ++s) {
    for (int l = 0; l < nL; ++l) {
      double sep, rej;
      metrics(hE[s][l], hP[s][l], sep, rej);
      printf("%-6s %-6d %6.0f %6.0f %8.2f %8.2f %8.3f %8.2f\n", sliceName[s], l,
             hE[s][l]->GetEntries(), hP[s][l]->GetEntries(),
             hE[s][l]->GetMean(), hP[s][l]->GetMean(), sep, rej);
    }
  }

  // ----- per-track table -----
  printf("\nPER TRACK (mean over valid layers, >= 4 layers)\n");
  printf("%-10s %6s %6s %8s %8s %8s %8s\n", "slice", "N_e", "N_pi", "mean_e", "mean_pi", "sep", "rej@90");
  for (int s = 0; s <= nS; ++s) {
    double sep, rej;
    metrics(hAvgE[s], hAvgP[s], sep, rej);
    printf("%-10s %6.0f %6.0f %8.2f %8.2f %8.3f %8.2f\n", sliceName[s],
           hAvgE[s]->GetEntries(), hAvgP[s]->GetEntries(),
           hAvgE[s]->GetMean(), hAvgP[s]->GetMean(), sep, rej);
  }
  printf("\nDesign goal for reference: rejection ~100 at 90%% electron efficiency with 6 layers.\n");
  printf("With ~200 electrons, treat every number above as having a large statistical error.\n\n");

  // ----- plots -----
  auto prep = [](TH1D* e, TH1D* p) {
    if (e->Integral() > 0) e->Scale(1.0 / e->Integral());
    if (p->Integral() > 0) p->Scale(1.0 / p->Integral());
    e->SetLineColor(kRed + 1);  e->SetLineWidth(2);
    p->SetLineColor(kBlue + 1); p->SetLineWidth(2);
    double m = std::max(e->GetMaximum(), p->GetMaximum());
    e->SetMaximum(1.25 * m);
  };

  // one canvas per slice: 6 layer pads
  for (int s = 0; s < nS; ++s) {
    TCanvas* c = new TCanvas(Form("c_slice%d", s), Form("%s per layer", sliceName[s]), 1500, 900);
    c->Divide(3, 2);
    for (int l = 0; l < nL; ++l) {
      c->cd(l + 1);
      prep(hE[s][l], hP[s][l]);
      hE[s][l]->Draw("hist");
      hP[s][l]->Draw("hist same");
      if (l == 0) {
        TLegend* lg = new TLegend(0.55, 0.72, 0.88, 0.88);
        lg->AddEntry(hE[s][l], "electrons (gamma)", "l");
        lg->AddEntry(hP[s][l], "pions (K0s)", "l");
        lg->Draw();
      }
    }
    c->SaveAs(Form("/home/madhvafakare/alice/my-analysis/trd/plots/trd_%s_per_layer.png", sliceName[s]));
  }

  // per-track averages, 4 pads (Q0, Q1, Q2, total)
  TCanvas* ct = new TCanvas("c_track", "per track", 1200, 900);
  ct->Divide(2, 2);
  for (int s = 0; s <= nS; ++s) {
    ct->cd(s + 1);
    prep(hAvgE[s], hAvgP[s]);
    hAvgE[s]->Draw("hist");
    hAvgP[s]->Draw("hist same");
    if (s == 0) {
      TLegend* lg = new TLegend(0.55, 0.72, 0.88, 0.88);
      lg->AddEntry(hAvgE[s], "electrons", "l");
      lg->AddEntry(hAvgP[s], "pions", "l");
      lg->Draw();
    }
  }
  ct->SaveAs("/home/madhvafakare/alice/my-analysis/trd/plots/trd_per_track_average.png");

  printf("Plots saved in ./plots/\n");
}