#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>

// Step 2c: P(Q|e) and P(Q|pi) for every window and every layer,
//          fQ00 ... fQ25, in 1 < pT < 2 GeV/c.
// Run:    root -l lq_step2c.C
// Output: ../plots/step2c_all.png  (rows = Q0, Q1, Q2; columns = layer 0 to 5)

void lq_step2c() {

    // =========================================================================
    // 1. Open the data
    // =========================================================================
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) {
        return;
    }

    // =========================================================================
    // 2. Settings
    // =========================================================================
    const int   nW    = 3;      // windows Q0, Q1, Q2
    const int   nBins = 15;     // only ~50-80 electrons per histogram, so keep it coarse
    const float ptMin = 1.0;
    const float ptMax = 2.0;

    // =========================================================================
    // 3. Axis range per histogram
    //    For each branch: find the charge that 99% of the tracks are below,
    //    and end the axis 10% above that. (Using the maximum instead would let
    //    a few outliers squeeze everything into the first bins.)
    // =========================================================================
    float qMax[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {

            // which branch, and which tracks
            TString branch    = Form("fQ%d%d", w, l);
            TString selection = Form("%s>0 && fSpecies>=0 && fPt>=%f && fPt<=%f",
                                     branch.Data(), ptMin, ptMax);

            // fill a temporary histogram with that branch (no window)
            c->Draw(branch + ">>hTmp(400,0,400)", selection, "goff");
            TH1* hTmp = (TH1*)gDirectory->Get("hTmp");

            // ask: below which value are 99% of the entries?
            double p = 0.99;     // the question: 99%
            double q = 100.0;    // the answer box (100 is only a fallback)
            if (hTmp && hTmp->GetEntries() > 0) {
                hTmp->GetQuantiles(1, &q, &p);
            }

            qMax[w][l] = 1.1 * q;
        }
    }

    // =========================================================================
    // 4. Book the histograms: one electron and one pion histogram per
    //    window and layer
    // =========================================================================
    TH1F* hEl[nW][kNLayers];
    TH1F* hPi[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            TString title = Form("Q_{%d}, layer %d;Q_{%d} (a.u.);normalized counts", w, l, w);

            hEl[w][l] = new TH1F(Form("hEl_%d%d", w, l), title, nBins, 0.0, qMax[w][l]);
            hPi[w][l] = new TH1F(Form("hPi_%d%d", w, l), "",    nBins, 0.0, qMax[w][l]);
        }
    }

    // =========================================================================
    // 5. Loop over the tracks and fill
    // =========================================================================
    for (Long64_t i = 0; i < c->GetEntries(); i++) {
        c->GetEntry(i);

        // momentum window
        if (t.pt < ptMin || t.pt > ptMax) {
            continue;
        }

        for (int w = 0; w < nW; w++) {
            for (int l = 0; l < kNLayers; l++) {

                float q = t.q[w][l];

                // this window is missing (-1)
                if (q <= 0) {
                    continue;
                }

                // a charge above the axis goes into the last bin, not lost
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

    // =========================================================================
    // 6. Print how many electrons and pions are in each histogram
    // =========================================================================
    printf("\nentries per window and layer\n");

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            int nEl = (int)hEl[w][l]->GetEntries();
            int nPi = (int)hPi[w][l]->GetEntries();
            printf("Q%d layer %d:  %d electrons,  %d pions\n", w, l, nEl, nPi);
        }
    }

    // =========================================================================
    // 7. Normalise to area 1: counts become probabilities P(Q|e) and P(Q|pi)
    // =========================================================================
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

    // =========================================================================
    // 8. Draw: 6 columns (layers) x 3 rows (windows)
    // =========================================================================
    gStyle->SetOptStat(0);

    TCanvas* cv = new TCanvas("cv", "", 2400, 1200);
    cv->Divide(kNLayers, nW);

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hEl[0][0], "electrons", "l");
    leg->AddEntry(hPi[0][0], "pions",     "l");

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {

            // pads are numbered 1 to 18, row by row
            int pad = w * kNLayers + l + 1;
            cv->cd(pad);

            // leave room above the higher of the two curves
            double top = std::max(hEl[w][l]->GetMaximum(), hPi[w][l]->GetMaximum());
            hEl[w][l]->SetMaximum(1.2 * top);

            hEl[w][l]->SetLineColor(kBlue + 1);
            hEl[w][l]->SetLineWidth(2);
            hPi[w][l]->SetLineColor(kRed + 1);
            hPi[w][l]->SetLineWidth(2);

            hEl[w][l]->Draw("hist");
            hPi[w][l]->Draw("hist same");

            // legend only in the first pad
            if (pad == 1) {
                leg->Draw();
            }
        }
    }

    cv->SaveAs("../plots/step2c_all.png");
}