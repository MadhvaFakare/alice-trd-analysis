#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>

/*
   Step 4c: same as step 4b, plus the clamp.
            Every value read from a table is kept between kClamp and 1 - kClamp,
            so no single table can say "certainly a pion" (L = 0) or
            "certainly an electron" (L = 1). One empty bin then costs a track
            a little, instead of killing it.
            The combined likelihood is computed twice per track, without and
            with the clamp, so the effect is visible side by side.
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

/* ===== Helper: normalise a histogram to area 1 ===== */
void Normalise(TH1F* h) {
    double area = h->Integral();
    if (area > 0) {
        h->Scale(1.0 / area);
    }
}

void lq_step4c() {

    /* ===== 1. Open the data ===== */
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) {
        return;
    }

    /* ===== 2. Settings ===== */
    const int    nW        = 3;
    const int    nBins     = 15;
    const float  ptMin     = 1.0;
    const float  ptMax     = 2.0;
    const int    minLayers = 4;
    const double kClamp    = 0.02;   /* NEW IN 4c: a table may say at most 98% / at least 2% */

    /* ===== 3. Axis range per histogram (same as step 4b) ===== */
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

    /* ===== 4. Book the histograms (same as step 4b) ===== */
    TH1F* hEl[nW][kNLayers];
    TH1F* hPi[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            hEl[w][l] = new TH1F(Form("hEl_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
            hPi[w][l] = new TH1F(Form("hPi_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
        }
    }

    /* ===== 5. Fill the references from the training half (same as step 4b) ===== */
    for (Long64_t i = 0; i < c->GetEntries(); i++) {
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

    /* ===== 6. Normalise to area 1 (same as step 4b) ===== */
    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            Normalise(hEl[w][l]);
            Normalise(hPi[w][l]);
        }
    }

    /* ===== 7. Lookup tables from the smoothed references (same as step 4b) ===== */
    TH1F* lut[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            TH1F* sEl = SmoothCopy(hEl[w][l], Form("sEl_%d%d", w, l));
            TH1F* sPi = SmoothCopy(hPi[w][l], Form("sPi_%d%d", w, l));
            lut[w][l] = BuildLut(sEl, sPi, Form("lut_%d%d", w, l));
        }
    }

    /* ===== 8. Combine the test half, without and with the clamp ===== */
    TH1F* hRawEl   = new TH1F("hRawEl",   "without clamp;L_{e};normalized counts",                     20, 0.0, 1.0);
    TH1F* hRawPi   = new TH1F("hRawPi",   "",                                                           20, 0.0, 1.0);
    TH1F* hClampEl = new TH1F("hClampEl", Form("with clamp %.2f;L_{e};normalized counts", kClamp),     20, 0.0, 1.0);
    TH1F* hClampPi = new TH1F("hClampPi", "",                                                           20, 0.0, 1.0);

    int nEl            = 0;
    int nPi            = 0;
    int nElKilledRaw   = 0;   /* electrons with L_e < 0.01, without clamp */
    int nElKilledClamp = 0;   /* electrons with L_e < 0.01, with clamp */
    int nUndefined     = 0;

    for (Long64_t i = 0; i < c->GetEntries(); i++) {
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

        double sumRaw   = 0.0;
        double sumClamp = 0.0;

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

                /* without clamp, as in step 4b */
                sumRaw += std::log((1.0 - L) / L);

                /* NEW IN 4c: keep L between kClamp and 1 - kClamp */
                double Lc = std::max(kClamp, std::min(1.0 - kClamp, L));
                sumClamp += std::log((1.0 - Lc) / Lc);
            }
        }

        if (std::isnan(sumRaw)) {
            nUndefined++;
            continue;
        }

        double LeRaw   = 1.0 / (1.0 + std::exp(sumRaw));
        double LeClamp = 1.0 / (1.0 + std::exp(sumClamp));
        if (LeRaw >= 1.0) {
            LeRaw = 0.999999;
        }
        if (LeClamp >= 1.0) {
            LeClamp = 0.999999;
        }

        if (t.species == kSpeciesElectron) {
            hRawEl->Fill(LeRaw);
            hClampEl->Fill(LeClamp);
            nEl++;
            if (LeRaw < 0.01) {
                nElKilledRaw++;
            }
            if (LeClamp < 0.01) {
                nElKilledClamp++;
            }
        } else {
            hRawPi->Fill(LeRaw);
            hClampPi->Fill(LeClamp);
            nPi++;
        }
    }

    /* ===== 9. Print ===== */
    printf("\ntest half: %d electrons, %d pions\n", nEl, nPi);
    printf("electrons with L_e < 0.01, without clamp: %d  (%.0f%%)\n",
           nElKilledRaw, nEl > 0 ? 100.0 * nElKilledRaw / nEl : 0.0);
    printf("electrons with L_e < 0.01, with clamp:    %d  (%.0f%%)\n",
           nElKilledClamp, nEl > 0 ? 100.0 * nElKilledClamp / nEl : 0.0);
    printf("tracks skipped as undefined (no clamp):   %d\n", nUndefined);

    /* ===== 10. Draw: left without clamp, right with clamp ===== */
    Normalise(hRawEl);
    Normalise(hRawPi);
    Normalise(hClampEl);
    Normalise(hClampPi);

    gStyle->SetOptStat(0);

    TCanvas* cv = new TCanvas("cv", "", 1600, 700);
    cv->Divide(2, 1);

    TH1F* elec[2] = {hRawEl, hClampEl};
    TH1F* pion[2] = {hRawPi, hClampPi};

    TLegend* leg = new TLegend(0.40, 0.78, 0.70, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hRawEl, "electrons", "l");
    leg->AddEntry(hRawPi, "pions",     "l");

    for (int k = 0; k < 2; k++) {
        cv->cd(k + 1)->SetLogy();

        elec[k]->SetMinimum(1e-3);
        elec[k]->SetMaximum(2.0);
        elec[k]->SetLineColor(kBlue + 1);
        elec[k]->SetLineWidth(2);
        pion[k]->SetLineColor(kRed + 1);
        pion[k]->SetLineWidth(2);

        elec[k]->Draw("hist");
        pion[k]->Draw("hist same");
        leg->Draw();
    }

    cv->SaveAs("../plots/step4c_clamp.png");
}