#include "trd_io.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void lq_step2(){
    /* Open */
    Track t;
    TChain* c = OpenTrdTree(t);
    if (!c) return;
    /* Settings */
    const int layer = 0; // which layer to look at from 0 to 5
    const int nBins = 20; 
    const float qMax = 200.0; // 99% of all totals are below this range
    /* Histogram made before the loop */
    TH1F* hEl = new TH1F("hEl","layer 0, 1 < p_{T} < 2 GeV/c; Q_{0}+Q_{2}+Q_{2} (a.u.); normalized counts", nBins, 0.0, qMax);
    TH1F* hPi = new TH1F("hPi","",nBins,0.0,qMax);
    /* Loop */
    for (Long64_t i = 0; i < c->GetEntries(); i++){
        c->GetEntry(i);
        if (t.pt < 1.0 || t.pt > 2.0) continue; // momentum window
        float q0 = t.q[0][layer];
        float q1 = t.q[1][layer];
        float q2 = t.q[2][layer];
        if (q0 <= 0 || q1 <= 0 || q2 <= 0) continue; // a window is missing (-1)
        float qTot = q0+q1+q2;
        if (t.species == kSpeciesElectron) hEl->Fill(qTot);
        if (t.species == kSpeciesPion) hPi->Fill(qTot);
    }
    /* Result */
    printf("electrons in layer %d:  %d\n",layer, (int)hEl->GetEntries());
    printf("pions in layer %d:  %d\n",layer, (int)hPi->GetEntries());
    /* Normalize to area 1: counts become probabilities P(Q|e) and P(Q|pi) */
    if (hEl->Integral() > 0) hEl->Scale(1.0/hEl->Integral());
    if (hPi->Integral() > 0) hPi->Scale(1.0/hPi->Integral());
    /* Draw both on same canvas */
    gStyle->SetOptStat(0);
    TCanvas* cv = new TCanvas("cv","",800,600);
    double top = std::max(hEl->GetMaximum(), hPi->GetMaximum());
    hEl->SetMaximum(1.2*top);

    hEl->SetLineColor(kBlue + 1); hEl->SetLineWidth(2); 
    hPi->SetLineColor(kRed + 1);  hPi->SetLineWidth(2);
    hEl->Draw("hist"); hPi->Draw("hist same");
    
    TLegend* leg = new TLegend(0.60,0.75,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hEl, "electrons", "l");
    leg->AddEntry(hPi, "pions", "l");
    leg->Draw();

    cv->SaveAs("../plots/step2_layer0.png");

}