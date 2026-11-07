// TRD Q0 charge: V0-tagged electrons vs pions
//
// Reads the PIDStudy task output histograms and quantifies the separation
// instead of only overlaying the two distributions.
//
// Usage:  root -l 'trd_q0_comparison2.C'
//         root -l 'trd_q0_comparison2.C("/path/to/AnalysisResults.root")'

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TMath.h"

void trd_q0_comparison2(
    const char *inFile = "/home/madhvafakare/alice/PbPb2023pass5/trd/AnalysisResults.root",
    const char *outPng = "/home/madhvafakare/alice/my-analysis/plots/trd_q0_comparison2.png") {

    TFile *f = TFile::Open(inFile);
    if (!f || f->IsZombie()) {
        printf("Error opening file: %s\n", inFile);
        return;
    }

    // NOTE: these are histograms inside the task's output directory,
    // not branches of a TTree. No flat tree is written yet.
    TH1F *hEl0 = (TH1F*)f->Get("t-r-d-p-i-d-study/Q0_El");
    TH1F *hPi0 = (TH1F*)f->Get("t-r-d-p-i-d-study/Q0_Pi");
    if (!hEl0 || !hPi0) {
        printf("Could not find Q0_El / Q0_Pi in t-r-d-p-i-d-study/\n");
        f->ls();
        return;
    }

    // ---- capture the numbers BEFORE any scaling ----------------------------
    const double nEl    = hEl0->GetEntries();
    const double nPi    = hPi0->GetEntries();
    const double meanEl = hEl0->GetMean();
    const double meanPi = hPi0->GetMean();
    const double rmsEl  = hEl0->GetRMS();
    const double rmsPi  = hPi0->GetRMS();

    printf("\n--- TRD Q0, V0-tagged samples ---\n");
    printf("Electrons: N = %.0f, mean = %.2f, RMS = %.2f\n", nEl, meanEl, rmsEl);
    printf("Pions    : N = %.0f, mean = %.2f, RMS = %.2f\n", nPi, meanPi, rmsPi);

    // Separation power: difference of means in units of the average width.
    const double sep = (rmsEl + rmsPi) > 0
                     ? (meanEl - meanPi) / (0.5 * (rmsEl + rmsPi)) : 0.0;
    printf("Separation power (dMean / <RMS>) = %.3f\n", sep);

    // ---- electron efficiency vs pion contamination -------------------------
    // Cut is "keep tracks with Q0 above threshold". Scan bin edges from the
    // top down and report the pion efficiency at fixed electron efficiency.
    const int nb = hEl0->GetNbinsX();
    const double totEl = hEl0->Integral(1, nb);
    const double totPi = hPi0->Integral(1, nb);
    if (totEl <= 0 || totPi <= 0) {
        printf("One of the histograms is empty, skipping the efficiency scan.\n");
    } else {
        const double targets[] = {0.95, 0.90, 0.80, 0.50};
        printf("\n%-12s %10s %12s %12s %12s\n",
               "target eff_e", "Q0 cut", "eff_e", "eff_pi", "rejection");
        for (double target : targets) {
            int bcut = -1;
            for (int b = nb; b >= 1; --b) {
                if (hEl0->Integral(b, nb) / totEl >= target) { bcut = b; break; }
            }
            if (bcut < 1) continue;
            const double cut   = hEl0->GetXaxis()->GetBinLowEdge(bcut);
            const double effEl = hEl0->Integral(bcut, nb) / totEl;
            const double effPi = hPi0->Integral(bcut, nb) / totPi;
            printf("%-12.2f %10.2f %12.4f %12.4f %12s\n", target, cut, effEl, effPi,
                   effPi > 0 ? Form("%.2f", 1.0 / effPi) : "inf");
        }
        printf("\nStatistical caveat: the electron sample has only %.0f entries,\n", nEl);
        printf("so these numbers carry a relative uncertainty of order %.1f%%.\n",
               100.0 / TMath::Sqrt(TMath::Max(nEl, 1.0)));
    }

    // ---- normalised overlay, linear and log --------------------------------
    TH1F *hEl = (TH1F*)hEl0->Clone("hEl_norm");
    TH1F *hPi = (TH1F*)hPi0->Clone("hPi_norm");
    hEl->Sumw2();
    hPi->Sumw2();
    if (hEl->Integral() > 0) hEl->Scale(1.0 / hEl->Integral());
    if (hPi->Integral() > 0) hPi->Scale(1.0 / hPi->Integral());

    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("cQ0Comparison", "TRD Q0: electrons vs pions", 1400, 600);
    c->Divide(2, 1);

    hEl->SetLineColor(kRed);   hEl->SetLineWidth(2);
    hPi->SetLineColor(kBlue);  hPi->SetLineWidth(2);
    hEl->SetTitle("TRD Q0 charge: electrons vs pions (V0-tagged);Q0 (a.u.);Normalised entries");

    for (int pad = 1; pad <= 2; ++pad) {
        c->cd(pad);
        if (pad == 2) gPad->SetLogy();
        hEl->Draw("HIST");
        hPi->Draw("HIST SAME");

        TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
        leg->SetBorderSize(1);
        leg->SetFillColor(kWhite);
        TLegendEntry *e1 = leg->AddEntry(hEl,
            Form("Electrons (mean = %.2f, N = %.0f)", meanEl, nEl), "l");
        e1->SetTextColor(kRed);
        TLegendEntry *e2 = leg->AddEntry(hPi,
            Form("Pions (mean = %.2f, N = %.0f)", meanPi, nPi), "l");
        e2->SetTextColor(kBlue);
        leg->Draw();
    }

    c->SaveAs(outPng);
    printf("\nWritten: %s\n", outPng);
}