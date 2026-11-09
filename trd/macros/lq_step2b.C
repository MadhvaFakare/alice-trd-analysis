#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>

void lq_step2b(){
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) return;

/* Settings */
const int layer = 0; // which layer to look at from 0 to 5
const int nBins = 20;
/* Four charges per layer: Q0, Q1, Q2 and their sum */
const int nQ = 4;
const char* qName[nQ] = {"Q_{0}", "Q_{1}", "Q_{2}", "Q_{0}+Q_{1}+Q_{2}"};
const float qMax[nQ] = {160.0, 180.0, 85.0, 200.0}; // the maxima from the data

/* Histogram made before loops */
TH1F* hEl[nQ];
TH1F* hPi[nQ];
for (int k = 0; k < nQ; k++){
    hEl[k] = new TH1F(Form("hEl%d", k ), Form("layer %d, 1 < p_{T} < 2 GeV/c; %s(a.u.)",layer,qName[k]), nBins, 0.0, qMax[k]);
    hPi[k] = new TH1F(Form("hPi%d",k), "",nBins,0.0,qMax[k]);
}

/* Loop */

for (Long64_t i = 0; i < c->GetEntries(); i++){
    c->GetEntry(i);
    if (t.pt < 1.0 || t.pt > 2.0) continue; // Momentum Window
    float q[nQ];
    q[0] = t.q[0][layer];
    q[1] = t.q[1][layer];
    q[2] = t.q[2][layer];
    if (q[0]<= 0 || q[1]<= 0 || q[2]<= 0) continue; // a window is missing (-1)
    q[3] = q[0] + q[1] + q[2];

    for (int k = 0; k < nQ; k++){
        if (t.species == kSpeciesElectron) hEl[k]->Fill(q[k]);
        if (t.species == kSpeciesPion) hPi[k]->Fill(q[k]);
    }
}
/* Results every Histogram holds the same track. One is enough */
printf("Electron in layer %d: %d\n",layer, (int)hEl[3]->GetEntries());
printf("Pions in layer %d: %d\n",layer, (int)hPi[3]->GetEntries());

/* Normalize to Area 1: counts becomes probabilities P(Q|e) and P(Q|pi) */
for (int k = 0; k < nQ; k++){
    if(hEl[k]->Integral() > 0) hEl[k]->Scale(1.0/hEl[k]->Integral());
    if(hPi[k]->Integral() > 0) hPi[k]->Scale(1.0/hPi[k]->Integral());
}

/* Draw : One pad per charge 2x2 canvas */

gStyle->SetOptStat(0);
TCanvas* cv = new TCanvas("cv", "", 1200,900);
cv->Divide(2,2);

TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
leg->SetBorderSize(0);
leg->AddEntry(hEl[0], "electrons", "l");
leg->AddEntry(hPi[0], "pions", "l");

for (int k = 0; k < nQ; k++){
    cv->cd(k+1);
    double top = std::max(hEl[k]->GetMaximum(), hPi[k]->GetMaximum());
    hEl[k]->SetMaximum(1.2*top);

    hEl[k]->SetLineColor(kBlue +1); hEl[k]->SetLineWidth(2);
    hPi[k]->SetLineColor(kRed + 1); hPi[k]->SetLineWidth(2);
    hEl[k]->Draw("hist");
    hPi[k]->Draw("hist same");
    if (k == 0) leg->Draw();
}
cv->SaveAs(Form("../plots/step2b_layer%d.png", layer));
}