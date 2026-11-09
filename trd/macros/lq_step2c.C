#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>

// Step 2c: P(Q|e) and P(Q|pi) for every window and every layer,
//          fQ00 ... fQ25, in 1 < pT < 2 GeV/c.
// Run:  root -l lq_step2c.C
// Output: ../plots/step2c_all.png  (rows = Q0, Q1, Q2; columns = layer 0 to 5)

void lq_step2c(){
    /* Open */
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) return;

    /* Settings */
    const int   nW    = 3;       // windows Q0, Q1, Q2
    const int   nBins = 15;      // only ~50-80 electrons per histogram, so keep it coarse
    const float ptMin = 1.0;
    const float ptMax = 2.0;

    /* Range per histogram: 1.1 x the 99% point of that branch, for the tagged
       tracks in the pT window. Using the maximum instead squeezes everything
       into the first bins, because a few outliers set the maximum. */
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

    /* Histograms made before the loop: one pair per window and layer */
    TH1F* hEl[nW][kNLayers];
    TH1F* hPi[nW][kNLayers];
    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            hEl[w][l] = new TH1F(Form("hEl_%d%d", w, l),
                                 Form("Q_{%d}, layer %d;Q_{%d} (a.u.);normalized counts", w, l, w),
                                 nBins, 0.0, qMax[w][l]);
            hPi[w][l] = new TH1F(Form("hPi_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
        }
    }

    /* Loop */
    for (Long64_t i = 0; i < c->GetEntries(); i++){
        c->GetEntry(i);
        if (t.pt < ptMin || t.pt > ptMax) continue;   // momentum window

        for (int w = 0; w < nW; w++){
            for (int l = 0; l < kNLayers; l++){
                float q = t.q[w][l];
                if (q <= 0) continue;                              // this window is missing (-1)
                if (q >= qMax[w][l]) q = 0.999 * qMax[w][l];       // tail goes into the last bin
                if (t.species == kSpeciesElectron) hEl[w][l]->Fill(q);
                if (t.species == kSpeciesPion)     hPi[w][l]->Fill(q);
            }
        }
    }

    /* Result: counts per histogram, electrons (pions) */
    printf("\nentries per window and layer: electrons (pions)\n");
    printf("     ");
    for (int l = 0; l < kNLayers; l++) printf("   layer %d   ", l);
    printf("\n");
    for (int w = 0; w < nW; w++){
        printf("Q%d   ", w);
        for (int l = 0; l < kNLayers; l++)
            printf("%4d (%5d) ", (int)hEl[w][l]->GetEntries(), (int)hPi[w][l]->GetEntries());
        printf("\n");
    }

    /* Normalize to area 1: counts become probabilities P(Q|e) and P(Q|pi) */
    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            if (hEl[w][l]->Integral() > 0) hEl[w][l]->Scale(1.0/hEl[w][l]->Integral());
            if (hPi[w][l]->Integral() > 0) hPi[w][l]->Scale(1.0/hPi[w][l]->Integral());
        }
    }

    /* Draw: 6 columns (layers) x 3 rows (windows) */
    gStyle->SetOptStat(0);
    TCanvas* cv = new TCanvas("cv", "", 2400, 1200);
    cv->Divide(kNLayers, nW);

    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hEl[0][0], "electrons", "l");
    leg->AddEntry(hPi[0][0], "pions", "l");

    for (int w = 0; w < nW; w++){
        for (int l = 0; l < kNLayers; l++){
            cv->cd(w * kNLayers + l + 1);                  // pads are numbered 1 to 18, row by row
            double top = std::max(hEl[w][l]->GetMaximum(), hPi[w][l]->GetMaximum());
            hEl[w][l]->SetMaximum(1.2 * top);

            hEl[w][l]->SetLineColor(kBlue + 1); hEl[w][l]->SetLineWidth(2);
            hPi[w][l]->SetLineColor(kRed + 1);  hPi[w][l]->SetLineWidth(2);
            hEl[w][l]->Draw("hist"); hPi[w][l]->Draw("hist same");
            if (w == 0 && l == 0) leg->Draw();
        }
    }

    cv->SaveAs("../plots/step2c_all.png");
}