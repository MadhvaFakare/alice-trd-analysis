#include "trd_io.h"
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>

/*
   Step 5: the efficiency curve and the pion rejection factor (Sec. 5.3.1).
           Same as step 4c with the clamp, then:
             - move a threshold down from L_e = 1, bin by bin
             - electron efficiency = fraction of electrons above the threshold
             - pion efficiency     = fraction of pions above the threshold
             - stop when the electron efficiency reaches 90%: the working point
             - rejection factor = 1 / pion efficiency
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

/* ===== Helper: normalise a histogram to area 1 (same as step 4c) ===== */
void Normalise(TH1F* h) {
    double area = h->Integral();
    if (area > 0) {
        h->Scale(1.0 / area);
    }
}

void lq_step5() {

    /* ===== 1. Open the data ===== */
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) {
        return;
    }

    /* ===== 2. Settings ===== */
    const int    nW         = 3;
    const int    nBins      = 15;
    const float  ptMin      = 1.0;
    const float  ptMax      = 2.0;
    const int    minLayers  = 4;
    const double kClamp     = 0.02;
    const double targetEff  = 0.90;   /* NEW IN 5: electron efficiency of the working point */
    const int    nLeBins    = 200;    /* NEW IN 5: fine bins for L_e, so the scan has small steps */

    /* ===== 3. Axis range per histogram (same as step 4c) ===== */
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

    /* ===== 4. Book the histograms (same as step 4c) ===== */
    TH1F* hEl[nW][kNLayers];
    TH1F* hPi[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            hEl[w][l] = new TH1F(Form("hEl_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
            hPi[w][l] = new TH1F(Form("hPi_%d%d", w, l), "", nBins, 0.0, qMax[w][l]);
        }
    }

    /* ===== 5. Fill the references from the training half (same as step 4c) ===== */
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

    /* ===== 6. Normalise to area 1 (same as step 4c) ===== */
    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            Normalise(hEl[w][l]);
            Normalise(hPi[w][l]);
        }
    }

    /* ===== 7. Lookup tables from the smoothed references (same as step 4c) ===== */
    TH1F* lut[nW][kNLayers];

    for (int w = 0; w < nW; w++) {
        for (int l = 0; l < kNLayers; l++) {
            TH1F* sEl = SmoothCopy(hEl[w][l], Form("sEl_%d%d", w, l));
            TH1F* sPi = SmoothCopy(hPi[w][l], Form("sPi_%d%d", w, l));
            lut[w][l] = BuildLut(sEl, sPi, Form("lut_%d%d", w, l));
        }
    }

    /* ===== 8. Combine the test half, with the clamp (same as step 4c) ===== */
    TH1F* hLeEl = new TH1F("hLeEl", "", nLeBins, 0.0, 1.0);
    TH1F* hLePi = new TH1F("hLePi", "", nLeBins, 0.0, 1.0);

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

        double sum = 0.0;
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
                L = std::max(kClamp, std::min(1.0 - kClamp, L));
                sum += std::log((1.0 - L) / L);
            }
        }

        double Le = 1.0 / (1.0 + std::exp(sum));
        if (Le >= 1.0) {
            Le = 0.999999;
        }

        if (t.species == kSpeciesElectron) {
            hLeEl->Fill(Le);
        } else {
            hLePi->Fill(Le);
        }
    }

    /* ===== 9. NEW IN 5: the efficiency scan =====
       Start at the top bin (L_e close to 1) and move down one bin at a time.
       After each step, accE and accP are the fractions of electrons and pions
       above the current threshold. */
    double nE = hLeEl->Integral();
    double nP = hLePi->Integral();

    TGraph* curve = new TGraph();
    double accE = 0.0;
    double accP = 0.0;

    double thr   = -1.0;   /* working point: threshold */
    double epsE  = -1.0;   /* working point: electron efficiency */
    double epsPi = -1.0;   /* working point: pion efficiency */

    for (int b = nLeBins; b >= 1; b--) {
        accE += hLeEl->GetBinContent(b) / nE;
        accP += hLePi->GetBinContent(b) / nP;

        /* a log axis cannot show 0, so only points with some pions */
        if (accP > 0) {
            curve->SetPoint(curve->GetN(), accE, accP);
        }

        /* the first time 90% of the electrons are above: working point */
        if (epsE < 0 && accE >= targetEff) {
            thr   = hLeEl->GetBinLowEdge(b);
            epsE  = accE;
            epsPi = accP;
        }
    }

    /* ===== 10. Print ===== */
    double rejection = (epsPi > 0) ? 1.0 / epsPi : -1.0;

    /* statistical uncertainty on the pion efficiency (binomial),
       and what it means for the rejection */
    double errPi  = (nP > 0) ? std::sqrt(epsPi * (1.0 - epsPi) / nP) : 0.0;
    double errRej = (epsPi > 0) ? errPi / (epsPi * epsPi) : 0.0;

    printf("\ntest half: %.0f electrons, %.0f pions\n", nE, nP);
    printf("working point:  threshold L_e > %.3f\n", thr);
    printf("  electron efficiency  %.3f\n", epsE);
    printf("  pion efficiency      %.3f +- %.3f\n", epsPi, errPi);
    printf("  PION REJECTION       %.2f +- %.2f\n", rejection, errRej);
    printf("Felix, LQ3D, 1-2 GeV/c: 3.36     design goal: 100\n");
    printf("note: one electron is %.1f%% of the electron efficiency\n", 100.0 / nE);

    /* ===== 11. Draw: left the likelihoods with the threshold, right the curve ===== */
    gStyle->SetOptStat(0);

    TCanvas* cv = new TCanvas("cv", "", 1600, 700);
    cv->Divide(2, 1);

    /* left: L_e for electrons and pions, in 20 bins for display */
    cv->cd(1)->SetLogy();
    TH1F* dEl = (TH1F*)hLeEl->Clone("dEl");
    TH1F* dPi = (TH1F*)hLePi->Clone("dPi");
    dEl->Rebin(10);
    dPi->Rebin(10);
    Normalise(dEl);
    Normalise(dPi);

    dEl->SetTitle("combined likelihood, test half;L_{e};normalized counts");
    dEl->SetMinimum(1e-3);
    dEl->SetMaximum(2.0);
    dEl->SetLineColor(kBlue + 1);
    dEl->SetLineWidth(2);
    dPi->SetLineColor(kRed + 1);
    dPi->SetLineWidth(2);
    dEl->Draw("hist");
    dPi->Draw("hist same");

    TLine* cut = new TLine(thr, 1e-3, thr, 2.0);
    cut->SetLineStyle(2);
    cut->Draw();

    TLegend* leg1 = new TLegend(0.40, 0.78, 0.70, 0.88);
    leg1->SetBorderSize(0);
    leg1->AddEntry(dEl, "electrons", "l");
    leg1->AddEntry(dPi, "pions",     "l");
    leg1->Draw();

    /* right: pion efficiency against electron efficiency */
    cv->cd(2)->SetLogy();
    gPad->DrawFrame(0.0, 1e-3, 1.0, 1.2,
                    Form("LQ3D, %.0f < p_{T} < %.0f GeV/c;electron efficiency #epsilon_{e};pion efficiency #epsilon_{#pi}",
                         ptMin, ptMax));
    curve->SetLineColor(kBlue + 1);
    curve->SetLineWidth(2);
    curve->Draw("L same");

    TLine* target = new TLine(targetEff, 1e-3, targetEff, 1.2);
    target->SetLineStyle(2);
    target->Draw();

    TGraph* wp = new TGraph();
    wp->SetPoint(0, epsE, epsPi);
    wp->SetMarkerStyle(20);
    wp->SetMarkerColor(kRed + 1);
    wp->Draw("P same");

    TLegend* leg2 = new TLegend(0.15, 0.78, 0.60, 0.88);
    leg2->SetBorderSize(0);
    leg2->AddEntry(curve, "LQ3D", "l");
    leg2->AddEntry(wp, Form("#epsilon_{#pi}(90%%) = %.2f, rejection %.1f", epsPi, rejection), "p");
    leg2->Draw();

    cv->SaveAs("../plots/step5_efficiency.png");
}