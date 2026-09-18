// Figure 5.5 from the flat tree: Q0, Q1, Q2 in layer 0 for electrons and pions,
// in a momentum window. This is what the histograms could never give you,
// because it needs per-track momentum and per-layer charge together.
//
// Usage: root -l 'macros/trd_fig55.C'
//        root -l 'macros/trd_fig55.C(1.0, 2.0)'   // different momentum window

void trd_fig55(double pMin = 1.0, double pMax = 2.0,
               const char *inFile =
               "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root") {

    // ---- chain every DF_*/trdpid in the file ----
    TFile *f = TFile::Open(inFile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", inFile); return; }

    TChain *ch = new TChain("trdpid");
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;
        ch->Add(TString(inFile) + "/" + name + "/trdpid");
    }
    printf("Chained %lld entries\n", ch->GetEntries());

    // ---- the selections ----
    TString cutEl = Form("fSpecies==0 && fPt>%f && fPt<%f", pMin, pMax);
    TString cutPi = Form("fSpecies==1 && fPt>%f && fPt<%f", pMin, pMax);

    printf("Momentum window %.1f to %.1f GeV/c\n", pMin, pMax);
    printf("  electrons in window : %lld\n", ch->GetEntries(cutEl));
    printf("  pions in window     : %lld\n", ch->GetEntries(cutPi));

    const char *var[3]   = {"fQ00", "fQ10", "fQ20"};
    const char *label[3] = {"Q0", "Q1", "Q2"};

    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("cFig55", "Layer 0 charges", 1500, 500);
    c->Divide(3, 1);

    for (int i = 0; i < 3; ++i) {
        c->cd(i + 1);

        TH1F *hEl = new TH1F(Form("hEl%d", i), "", 60, 0, 120);
        TH1F *hPi = new TH1F(Form("hPi%d", i), "", 60, 0, 120);

        ch->Project(hEl->GetName(), var[i], cutEl);
        ch->Project(hPi->GetName(), var[i], cutPi);

        const double nEl = hEl->GetEntries(), nPi = hPi->GetEntries();
        const double mEl = hEl->GetMean(),    mPi = hPi->GetMean();
        printf("%s layer 0: electrons mean %.2f (N %.0f), pions mean %.2f (N %.0f)\n",
               label[i], mEl, nEl, mPi, nPi);

        if (hEl->Integral() > 0) hEl->Scale(1.0 / hEl->Integral());
        if (hPi->Integral() > 0) hPi->Scale(1.0 / hPi->Integral());

        hEl->SetLineColor(kBlue);  hEl->SetLineWidth(2);
        hPi->SetLineColor(kRed);   hPi->SetLineWidth(2);
        hEl->SetTitle(Form("%s layer 0, %.1f-%.1f GeV/c;%s (a.u.);normalised",
                           label[i], pMin, pMax, label[i]));

        // draw whichever is taller first so nothing is cut off
        if (hEl->GetMaximum() > hPi->GetMaximum()) {
            hEl->Draw("HIST"); hPi->Draw("HIST SAME");
        } else {
            hPi->Draw("HIST"); hEl->Draw("HIST SAME");
        }

        TLegend *leg = new TLegend(0.50, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hEl, Form("e^{#pm}  mean %.2f  (N %.0f)", mEl, nEl), "l");
        leg->AddEntry(hPi, Form("#pi^{#pm} mean %.2f  (N %.0f)", mPi, nPi), "l");
        leg->Draw();
    }

    c->SaveAs("/home/madhvafakare/alice/my-analysis/plots/trd_fig55.png");
    printf("\nWritten: plots/trd_fig55.png\n");
}