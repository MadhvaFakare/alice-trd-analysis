#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <algorithm>

// Step 3: one lookup table per window and layer (Eq. 5.4)
//             L(Q) = P(Q|e) / [ P(Q|e) + P(Q|pi) ]
//         built from the histograms of step 2c, in 1 < pT < 2 GeV/c.
// Run:  root -l lq_step3.C
// Output: ../plots/step3_luts.png  (rows = Q0, Q1, Q2; columns = layer 0 to 5)

void lq_step3(){
    /* Open */
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) return;

    /* Settings (same as step 2c) */
    const int   nW    = 3;
    const int   nBins = 15;
    const float ptMin = 1.0;
    const float ptMax = 2.0;

    /* Range per histogram: 1.1 x the 99% point (same as step 2c) */
    float qMax[nW][kNLayers];
    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            c->Draw(Form("fQ%d%d>>hTmp(400,0,400)", w, l),
                    Form("fQ%d%d>0 && fSpecies>=0 && fPt>=%f && fPt<=%f", w, l, ptMin, ptMax),
                    "goff");
            TH1* hTmp = (TH1*)gDirectory->Get("hTmp");
            double p = 0.99, q = 100.0;
            if (hTmp && hTmp->GetEntries() > 0) hTmp->GetQuantiles(1, &q, &p);
            qMax[w][l] = 1.1 * q;
        }
    }

    /* Histograms (same as step 2c) */
    TH1F* hEl[nW][kNLayers];
    TH1F* hPi[nW][kNLayers];
    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            hEl[w][l] = new TH1F(Form("hEl_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
            hPi[w][l] = new TH1F(Form("hPi_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
        }
    }

    /* Loop (same as step 2c) */
    for (Long64_t i = 0; i < c->GetEntries(); i++){
        c->GetEntry(i);
        if (t.pt < ptMin || t.pt > ptMax) continue;
        for (int w = 0; w < nW; w++){
            for (int l = 0; l < kNLayers; l++){
                float q = t.q[w][l];
                if (q <= 0) continue;
                if (q >= qMax[w][l]) q = 0.999 * qMax[w][l];
                if (t.species == kSpeciesElectron) hEl[w][l]->Fill(q);
                if (t.species == kSpeciesPion)     hPi[w][l]->Fill(q);
            }
        }
    }

    /* Normalize to area 1 (same as step 2c) */
    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            if (hEl[w][l]->Integral() > 0) hEl[w][l]->Scale(1.0/hEl[w][l]->Integral());
            if (hPi[w][l]->Integral() > 0) hPi[w][l]->Scale(1.0/hPi[w][l]->Integral());
        }
    }

    /* ===== NEW IN STEP 3: the lookup tables ===== */
    TH1F* hL[nW][kNLayers];
    int nZero[nW][kNLayers];          // bins where L = 0: pions present, no electrons
    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            hL[w][l] = new TH1F(Form("hL_%d%d", w, l),
                                Form("Q_{%d}, layer %d;Q_{%d} (a.u.);L", w, l, w),
                                nBins, 0.0, qMax[w][l]);
            nZero[w][l] = 0;
            for (int b = 1; b <= nBins; b++){
                double e = hEl[w][l]->GetBinContent(b);
                double p = hPi[w][l]->GetBinContent(b);
                double L = 0.5;                       // empty bin: no information
                if (e + p > 0) L = e / (e + p);
                if (e == 0 && p > 0) nZero[w][l]++;
                hL[w][l]->SetBinContent(b, L);
            }
        }
    }

    /* Result: how many bins per table say "certainly a pion" (L = 0) */
    printf("\nbins with L = 0 (pions but no electrons), out of %d\n", nBins);
    printf("     ");
    for (int l = 0; l < kNLayers; l++) printf(" layer %d", l);
    printf("\n");
    for (int w = 0; w < nW; w++){
        printf("Q%d   ", w);
        for (int l = 0; l < kNLayers; l++) printf("%8d", nZero[w][l]);
        printf("\n");
    }

    /* Draw: 6 columns (layers) x 3 rows (windows), with a line at L = 0.5 */
    gStyle->SetOptStat(0);
    TCanvas* cv = new TCanvas("cv", "", 2400, 1200);
    cv->Divide(kNLayers, nW);
    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            cv->cd(w * kNLayers + l + 1);
            hL[w][l]->SetMinimum(0.0);
            hL[w][l]->SetMaximum(1.0);
            hL[w][l]->SetLineColor(kBlack);
            hL[w][l]->SetLineWidth(2);
            hL[w][l]->Draw("hist");
            TLine* half = new TLine(0.0, 0.5, qMax[w][l], 0.5);
            half->SetLineStyle(2);
            half->Draw();
        }
    }

    cv->SaveAs("../plots/step3_luts.png");
}