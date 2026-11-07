// TRD charge slices Q0, Q1, Q2: which one separates electrons from pions?
//
// The TR photon is absorbed near the chamber entrance and appears at long
// drift times, so the three time slices should not discriminate equally.
// This runs the same efficiency scan on each slice and prints a comparison.
//
// Usage:  root -l 'trd_q_slices.C'

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPad.h"

// Returns pion efficiency at the requested electron efficiency, and fills
// the cut value actually used.
double scanAt(TH1F *hEl, TH1F *hPi, double target, double &cutOut, double &effElOut) {
    const int nb = hEl->GetNbinsX();
    const double totEl = hEl->Integral(1, nb);
    const double totPi = hPi->Integral(1, nb);
    if (totEl <= 0 || totPi <= 0) return -1;
    for (int b = nb; b >= 1; --b) {
        if (hEl->Integral(b, nb) / totEl >= target) {
            cutOut   = hEl->GetXaxis()->GetBinLowEdge(b);
            effElOut = hEl->Integral(b, nb) / totEl;
            return hPi->Integral(b, nb) / totPi;
        }
    }
    return -1;
}

void trd_q_slices(
    const char *inFile = "/home/madhvafakare/alice/PbPb2023pass5/trd/AnalysisResults.root",
    const char *outPng = "/home/madhvafakare/alice/my-analysis/trd/plots/trd_q_slices.png") {

    TFile *f = TFile::Open(inFile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", inFile); return; }

    const char *slice[3] = {"Q0", "Q1", "Q2"};
    TH1F *hEl[3], *hPi[3];

    for (int i = 0; i < 3; ++i) {
        hEl[i] = (TH1F*)f->Get(Form("t-r-d-p-i-d-study/%s_El", slice[i]));
        hPi[i] = (TH1F*)f->Get(Form("t-r-d-p-i-d-study/%s_Pi", slice[i]));
        if (!hEl[i] || !hPi[i]) { printf("Missing %s_El or %s_Pi\n", slice[i], slice[i]); return; }
    }

    printf("\n=== TRD charge slices: electron / pion separation ===\n");
    printf("%-6s %10s %10s %10s %10s %12s %10s %12s\n",
           "slice", "mean_e", "mean_pi", "RMS_e", "RMS_pi", "significance", "sep.power", "rej@90%");

    for (int i = 0; i < 3; ++i) {
        const double nE = hEl[i]->GetEntries(), nP = hPi[i]->GetEntries();
        const double mE = hEl[i]->GetMean(),    mP = hPi[i]->GetMean();
        const double rE = hEl[i]->GetRMS(),     rP = hPi[i]->GetRMS();

        // significance of the mean difference
        const double seE  = nE > 0 ? rE / TMath::Sqrt(nE) : 0;
        const double seP  = nP > 0 ? rP / TMath::Sqrt(nP) : 0;
        const double sig  = (seE + seP) > 0 ? (mE - mP) / TMath::Sqrt(seE*seE + seP*seP) : 0;
        const double sepP = (rE + rP) > 0 ? (mE - mP) / (0.5 * (rE + rP)) : 0;

        double cut = 0, effEl = 0;
        const double effPi = scanAt(hEl[i], hPi[i], 0.90, cut, effEl);
        const double rej   = (effPi > 0) ? 1.0 / effPi : -1;

        printf("%-6s %10.2f %10.2f %10.2f %10.2f %10.1f s %10.3f %12s\n",
               slice[i], mE, mP, rE, rP, sig, sepP,
               rej > 0 ? Form("%.2f", rej) : "n/a");
    }

    printf("\nAt 50%% electron efficiency:\n");
    printf("%-6s %10s %10s %12s\n", "slice", "Q cut", "eff_pi", "rejection");
    for (int i = 0; i < 3; ++i) {
        double cut = 0, effEl = 0;
        const double effPi = scanAt(hEl[i], hPi[i], 0.50, cut, effEl);
        printf("%-6s %10.2f %10.4f %12s\n", slice[i], cut, effPi,
               effPi > 0 ? Form("%.2f", 1.0/effPi) : "n/a");
    }

    printf("\nNote: QTot_El / QTot_Pi are not in the file, so the per-track sum\n");
    printf("of the three slices cannot be built from histograms alone.\n");
    printf("That combination needs the flat tree.\n");

    // ---- three-panel normalised overlay, log y ----------------------------
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("cSlices", "TRD charge slices", 1500, 500);
    c->Divide(3, 1);
    for (int i = 0; i < 3; ++i) {
        const double mE = hEl[i]->GetMean(), mP = hPi[i]->GetMean();
        const double nE = hEl[i]->GetEntries(), nP = hPi[i]->GetEntries();

        TH1F *e = (TH1F*)hEl[i]->Clone(Form("e%d", i));
        TH1F *p = (TH1F*)hPi[i]->Clone(Form("p%d", i));
        if (e->Integral() > 0) e->Scale(1.0 / e->Integral());
        if (p->Integral() > 0) p->Scale(1.0 / p->Integral());

        c->cd(i + 1);
        gPad->SetLogy();
        e->SetLineColor(kRed);  e->SetLineWidth(2);
        p->SetLineColor(kBlue); p->SetLineWidth(2);
        e->SetTitle(Form("%s: electrons vs pions (V0-tagged);%s (a.u.);Normalised entries",
                         slice[i], slice[i]));
        e->Draw("HIST");
        p->Draw("HIST SAME");

        TLegend *leg = new TLegend(0.50, 0.74, 0.88, 0.88);
        leg->SetBorderSize(1);
        leg->AddEntry(e, Form("e  (mean %.2f, N %.0f)", mE, nE), "l");
        leg->AddEntry(p, Form("#pi (mean %.2f, N %.0f)", mP, nP), "l");
        leg->Draw();
    }
    c->SaveAs(outPng);
    printf("\nWritten: %s\n", outPng);
}