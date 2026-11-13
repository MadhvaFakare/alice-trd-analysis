#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>

/*
   Step 4b: same as step 4, but honest.
            The tables are built from the even-numbered tracks (training half)
            and the combined likelihood is computed for the odd-numbered tracks
            (test half). No track is tested by a table it helped to build.
            Still no clamp: this step shows what the empty bins do to tracks
            the tables have never seen.
*/

/* ===== Helper: simple smoothing (same as step 3b) ===== */
TH1F* SmoothCopy(TH1F* h, const char* name) {
    TH1F* s = (TH1F*)h->Clone(name);
    int n = h->GetNbinsX();

    for (int b = 1; b <= n; b++) {
        double mid   = h->GetBinContent(b);
        double left  = (b > 1) ? h->GetBinContent(b - 1) : mid;
        double right = (b < n) ? h->GetBinContent(b + 1) : mid;
        s->SetBinContent(b, 0.25 * left + 0.5 * mid + 0.25 * right);
    }

    double area = s->Integral();
    if (area > 0) {
        s->Scale(1.0 / area);
    }
    return s;
}

/* ===== Helper: lookup table L = e / (e + pi) (same as step 3b) ===== */
TH1F* BuildLut(TH1F* hE, TH1F* hP, const char* name) {
    TH1F* lut = (TH1F*)hE->Clone(name);
    lut->Reset();

    for (int b = 1; b <= lut->GetNbinsX(); b++) {
        double e = hE->GetBinContent(b);
        double p = hP->GetBinContent(b);
        double L = 0.5;
        if (e + p > 0) {
            L = e / (e + p);
        }
        lut->SetBinContent(b, L);
    }
    return lut;
}

void lq_step4b() {

    /* ===== 1. Open the data ===== */
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) {
        return;
    }

    /* ===== 2. Settings ===== */
    const int   nW        = 3;
    const int   nBins     = 15;
    const float ptMin     = 1.0;
    const float ptMax     = 2.0;
    const int   minLayers = 4;

    /* ===== 3. Axis range per histogram (same as step 4) =====
       This uses all tagged tracks, but it only decides where the axis ends,
       not what the tables say, so it does not leak the answer to the test. */
    float qMax[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            TString branch    = Form("fQ%d%d", w, l);
            TString selection = Form("%s>0 && fSpecies>=0 && fPt>=%f && fPt<=%f",
                                     branch.Data(), ptMin, ptMax);

            c->Draw(branch + ">>hTmp(400,0,400)", selection, "goff");
            TH1* hTmp = (TH1*)gDirectory->Get("hTmp");

            double p = 0.99;
            double q = 100.0;
            if (hTmp && hTmp->GetEntries() > 0) {
                hTmp->GetQuantiles(1, &q, &p);
            }
            qMax[w][l] = 1.1 * q;
        }
    }

    /* ===== 4. Book the histograms (same as step 4) ===== */
    TH1F* hEl[nW][kNLayers];
    TH1F* hPi[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            hEl[w][l] = new TH1F(Form("hEl_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
            hPi[w][l] = new TH1F(Form("hPi_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
        }
    }

    /* ===== 5. Fill the references from the TRAINING half only ===== */
    for (Long64_t i = 0; i < c->GetEntries(); i++) {

        /* NEW IN 4b: odd-numbered tracks are kept for testing */
        if (i % 2 == 1) {
            continue;
        }

        c->GetEntry(i);

        if (t.pt < ptMin || t.pt > ptMax) {
            continue;
        }

        for (int w = 0; w < nW; w++) {
            for (int l = 0; l < kNLayers; l++) {
                float q = t.q[w][l];
                if (q <= 0) {
                    continue;
                }
                if (q >= qMax[w][l]) {
                    q = 0.999 * qMax[w][l];
                }
                if (t.species == kSpeciesElectron) {
                    hEl[w][l]->Fill(q);
                }
                if (t.species == kSpeciesPion) {
                    hPi[w][l]->Fill(q);
                }
            }
        }
    }

    /* ===== 6. Normalise to area 1 (same as step 4) ===== */
    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            double areaEl = hEl[w][l]->Integral();
            double areaPi = hPi[w][l]->Integral();
            if (areaEl > 0) {
                hEl[w][l]->Scale(1.0 / areaEl);
            }
            if (areaPi > 0) {
                hPi[w][l]->Scale(1.0 / areaPi);
            }
        }
    }

    /* ===== 7. Lookup tables from the smoothed references (same as step 4) ===== */
    TH1F* lut[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            TH1F* sEl = SmoothCopy(hEl[w][l], Form("sEl_%d%d", w, l));
            TH1F* sPi = SmoothCopy(hPi[w][l], Form("sPi_%d%d", w, l));
            lut[w][l] = BuildLut(sEl, sPi, Form("lut_%d%d", w, l));
        }
    }

    /* ===== 8. Combine the tables for every track of the TEST half ===== */
    TH1F* hLeEl = new TH1F("hLeEl", "combined likelihood, test half;L_{e};normalized counts", 20, 0.0, 1.0);
    TH1F* hLePi = new TH1F("hLePi", "", 20, 0.0, 1.0);

    int nEl         = 0;   /* electrons combined */
    int nPi         = 0;   /* pions combined */
    int nElWithZero = 0;   /* electrons that hit at least one L = 0 bin */
    int nElKilled   = 0;   /* electrons that end up with L_e < 0.01 */
    int nPiWithOne  = 0;   /* pions that hit at least one L = 1 bin (the mirror problem) */
    int nUndefined  = 0;   /* tracks that hit both L = 0 and L = 1: infinity minus infinity */

    for (Long64_t i = 0; i < c->GetEntries(); i++) {

        /* NEW IN 4b: only the test half */
        if (i % 2 == 0) {
            continue;
        }

        c->GetEntry(i);

        if (t.pt < ptMin || t.pt > ptMax) {
            continue;
        }
        if (t.species != kSpeciesElectron && t.species != kSpeciesPion) {
            continue;
        }

        int nLayers = 0;
        for (int l = 0; l < kNLayers; l++) {
            if (t.q[0][l] > 0 && t.q[1][l] > 0 && t.q[2][l] > 0) {
                nLayers++;
            }
        }
        if (nLayers < minLayers) {
            continue;
        }

        double sum     = 0.0;
        bool   hitZero = false;
        bool   hitOne  = false;

        for (int l = 0; l < kNLayers; l++) {
            if (t.q[0][l] <= 0 || t.q[1][l] <= 0 || t.q[2][l] <= 0) {
                continue;
            }
            for (int w = 0; w < nW; w++) {
                float q = t.q[w][l];
                if (q >= qMax[w][l]) {
                    q = 0.999 * qMax[w][l];
                }

                double L = lut[w][l]->GetBinContent(lut[w][l]->FindFixBin(q));
                if (L == 0.0) {
                    hitZero = true;
                }
                if (L == 1.0) {
                    hitOne = true;
                }
                sum += std::log((1.0 - L) / L);
            }
        }

        /* NEW IN 4b: with unseen tracks, one track can hit L = 0 in one table
           and L = 1 in another. Then the sum is +infinity - infinity, which is
           undefined (NaN). Count it and skip the track. */
        if (std::isnan(sum)) {
            nUndefined++;
            continue;
        }

        double Le = 1.0 / (1.0 + std::exp(sum));
        if (Le >= 1.0) {
            Le = 0.999999;
        }

        if (t.species == kSpeciesElectron) {
            hLeEl->Fill(Le);
            nEl++;
            if (hitZero) {
                nElWithZero++;
            }
            if (Le < 0.01) {
                nElKilled++;
            }
        } else {
            hLePi->Fill(Le);
            nPi++;
            if (hitOne) {
                nPiWithOne++;
            }
        }
    }

    /* ===== 9. Print ===== */
    printf("\ntest half: %d electrons, %d pions (at least %d complete layers)\n",
           nEl, nPi, minLayers);
    printf("electrons that hit an L = 0 bin in at least one table: %d\n", nElWithZero);
    printf("electrons with combined L_e < 0.01:                    %d  (%.0f%%)\n",
           nElKilled, nEl > 0 ? 100.0 * nElKilled / nEl : 0.0);
    printf("pions that hit an L = 1 bin in at least one table:     %d\n", nPiWithOne);
    printf("tracks skipped, L = 0 and L = 1 at once (undefined):    %d\n", nUndefined);

    /* ===== 10. Draw the combined likelihood for electrons and pions ===== */
    double areaEl = hLeEl->Integral();
    double areaPi = hLePi->Integral();
    if (areaEl > 0) {
        hLeEl->Scale(1.0 / areaEl);
    }
    if (areaPi > 0) {
        hLePi->Scale(1.0 / areaPi);
    }

    gStyle->SetOptStat(0);

    TCanvas* cv = new TCanvas("cv", "", 900, 700);
    cv->SetLogy();

    hLeEl->SetMinimum(1e-3);
    hLeEl->SetMaximum(2.0);
    hLeEl->SetLineColor(kBlue + 1);
    hLeEl->SetLineWidth(2);
    hLePi->SetLineColor(kRed + 1);
    hLePi->SetLineWidth(2);

    hLeEl->Draw("hist");
    hLePi->Draw("hist same");

    TLegend* leg = new TLegend(0.40, 0.78, 0.70, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hLeEl, "electrons", "l");
    leg->AddEntry(hLePi, "pions",     "l");
    leg->Draw();

    cv->SaveAs("../plots/step4b_test_half.png");
}