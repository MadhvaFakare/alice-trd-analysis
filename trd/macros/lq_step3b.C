#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <algorithm>

/*
   Step 3b: smooth the reference histograms before building the lookup tables.
            With only ~50-80 electrons per histogram, some bins are empty by
            chance, which gives L = 0 ("certainly a pion") or L = 1 ("certainly
            an electron") for no physical reason. Smoothing fills those holes.
            The plot shows each table before (grey) and after (black) smoothing.
*/

/* ===== Helper: simple smoothing =====
   Each bin becomes 1/4 of its left neighbour + 1/2 of itself + 1/4 of its
   right neighbour. A single empty bin between filled ones gets filled; a bin
   that had entries never becomes empty. At the edges, the missing neighbour
   is replaced by the bin itself. Afterwards the area is set back to 1. */
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

/* ===== Helper: lookup table L = e / (e + pi), bin by bin (Eq. 5.4) ===== */
TH1F* BuildLut(TH1F* hE, TH1F* hP, const char* name) {
    TH1F* lut = (TH1F*)hE->Clone(name);
    lut->Reset();

    for (int b = 1; b <= lut->GetNbinsX(); b++) {
        double e = hE->GetBinContent(b);
        double p = hP->GetBinContent(b);
        double L = 0.5;              /* empty bin: no information */
        if (e + p > 0) {
            L = e / (e + p);
        }
        lut->SetBinContent(b, L);
    }
    return lut;
}

/* ===== Helper: count bins with L = 0 or L = 1 ===== */
int CountBins(TH1F* lut, double value) {
    int n = 0;
    for (int b = 1; b <= lut->GetNbinsX(); b++) {
        if (lut->GetBinContent(b) == value) {
            n++;
        }
    }
    return n;
}

void lq_step3b() {

    /* ===== 1. Open the data ===== */
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) {
        return;
    }

    /* ===== 2. Settings (same as step 3) ===== */
    const int   nW    = 3;
    const int   nBins = 15;
    const float ptMin = 1.0;
    const float ptMax = 2.0;

    /* ===== 3. Axis range per histogram (same as step 3) ===== */
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

    /* ===== 4. Book the histograms (same as step 3) ===== */
    TH1F* hEl[nW][kNLayers];
    TH1F* hPi[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            hEl[w][l] = new TH1F(Form("hEl_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
            hPi[w][l] = new TH1F(Form("hPi_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
        }
    }

    /* ===== 5. Fill (same as step 3) ===== */
    for (Long64_t i = 0; i < c->GetEntries(); i++) {
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

    /* ===== 6. Normalise to area 1 (same as step 3) ===== */
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

    /* ===== 7. NEW: tables before and after smoothing ===== */
    TH1F* lutRaw[nW][kNLayers];
    TH1F* lutSmooth[nW][kNLayers];

    printf("\nbins with L = 0 and L = 1, out of %d: before -> after smoothing\n", nBins);

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {

            /* before: straight from the histograms, as in step 3 */
            lutRaw[w][l] = BuildLut(hEl[w][l], hPi[w][l], Form("lutRaw_%d%d", w, l));

            /* after: smooth both references first, then build the table */
            TH1F* sEl = SmoothCopy(hEl[w][l], Form("sEl_%d%d", w, l));
            TH1F* sPi = SmoothCopy(hPi[w][l], Form("sPi_%d%d", w, l));
            lutSmooth[w][l] = BuildLut(sEl, sPi, Form("lutSmooth_%d%d", w, l));

            int zeroBefore = CountBins(lutRaw[w][l],    0.0);
            int zeroAfter  = CountBins(lutSmooth[w][l], 0.0);
            int oneBefore  = CountBins(lutRaw[w][l],    1.0);
            int oneAfter   = CountBins(lutSmooth[w][l], 1.0);

            printf("Q%d layer %d:  L=0  %d -> %d    L=1  %d -> %d\n",
                   w, l, zeroBefore, zeroAfter, oneBefore, oneAfter);
        }
    }

    /* ===== 8. Draw: grey = before, black = after, dashed line at L = 0.5 ===== */
    gStyle->SetOptStat(0);

    TCanvas* cv = new TCanvas("cv", "", 2400, 1200);
    cv->Divide(kNLayers, nW);

    TLegend* leg = new TLegend(0.15, 0.75, 0.60, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(lutRaw[0][0],    "raw",      "l");
    leg->AddEntry(lutSmooth[0][0], "smoothed", "l");

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {

            int pad = w * kNLayers + l + 1;
            cv->cd(pad);

            lutRaw[w][l]->SetTitle(Form("Q_{%d}, layer %d;Q_{%d} (a.u.);L", w, l, w));
            lutRaw[w][l]->SetMinimum(0.0);
            lutRaw[w][l]->SetMaximum(1.0);
            lutRaw[w][l]->SetLineColor(kGray + 1);
            lutRaw[w][l]->SetLineWidth(1);

            lutSmooth[w][l]->SetLineColor(kBlack);
            lutSmooth[w][l]->SetLineWidth(2);

            lutRaw[w][l]->Draw("hist");
            lutSmooth[w][l]->Draw("hist same");

            TLine* half = new TLine(0.0, 0.5, qMax[w][l], 0.5);
            half->SetLineStyle(2);
            half->Draw();

            if (pad == 1) {
                leg->Draw();
            }
        }
    }

    cv->SaveAs("../plots/step3b_luts_smoothed.png");
}