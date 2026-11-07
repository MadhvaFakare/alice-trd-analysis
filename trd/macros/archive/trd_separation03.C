// trd_separation03.C
//
// Per-slice, per-layer TRD charge: electrons (red) vs pions (blue), 3x6 grid.
// Also prints how many tagged tracks fall in each pT window, and plots the
// pT distributions of electrons and pions on a second canvas.
//
// Charge value <= 0 means "no tracklet" (-1 marker) and is skipped.
// ptMin/ptMax only affect the 3x6 charge plots. Defaults (0, 1e9) mean "no cut".
//
// Usage:
//   root -l 'macros/trd_separation03.C'
//   root -l 'macros/trd_separation03.C("/path/trdpid_tree.root", 0.5, 1.0)'

void trd_separation03(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root",
                      double ptMin = 0.0, double ptMax = 1e9)
{
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) return;
    gStyle->SetOptStat(0);

    const int nSlice = 3, nLayer = 6;
    double xmax[nSlice] = {128, 128, 128};

    TH1D* hE[nSlice][nLayer];
    TH1D* hP[nSlice][nLayer];
    for (int s = 0; s < nSlice; s++) {
        for (int l = 0; l < nLayer; l++) {
            hE[s][l] = new TH1D(Form("hE_%d%d", s, l), Form("Q%d in layer %d;charge;a.u.", s, l), 64, 0, xmax[s]);
            hP[s][l] = new TH1D(Form("hP_%d%d", s, l), Form("Q%d in layer %d;charge;a.u.", s, l), 64, 0, xmax[s]);
        }
    }

    // pT windows for the counting table (last window is open-ended)
    const int nWin = 5;
    const double edge[nWin + 1] = {0.0, 0.5, 1.0, 2.0, 4.0, 1e9};
    long winE[nWin] = {0}, winP[nWin] = {0};

    TH1D* hPtE = new TH1D("hPtE", "pT of tagged tracks (layer 0 usable);pT (GeV/c);a.u.", 50, 0, 10);
    TH1D* hPtP = new TH1D("hPtP", "pT of tagged tracks (layer 0 usable);pT (GeV/c);a.u.", 50, 0, 10);

    long nTot = 0, nUnknown = 0, nE = 0, nP = 0;
    double ptSumE = 0, ptSumP = 0, ptLo = 1e9, ptHi = -1e9;

    TIter next(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;
        TDirectory* dir = (TDirectory*)key->ReadObj();
        TTree* t = (TTree*)dir->Get("trdpid");
        if (!t) continue;

        TLeaf* lSp = t->GetLeaf("fSpecies");
        TLeaf* lPt = t->GetLeaf("fPt");
        TLeaf* lQ[nSlice][nLayer];
        bool ok = (lSp != nullptr && lPt != nullptr);
        for (int s = 0; s < nSlice && ok; s++)
            for (int l = 0; l < nLayer; l++) {
                lQ[s][l] = t->GetLeaf(Form("fQ%d%d", s, l));
                if (!lQ[s][l]) { ok = false; break; }
            }
        if (!ok) continue;

        // only read the branches we need (much faster)
        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("fSpecies", 1);
        t->SetBranchStatus("fPt", 1);
        for (int s = 0; s < nSlice; s++)
            for (int l = 0; l < nLayer; l++)
                t->SetBranchStatus(Form("fQ%d%d", s, l), 1);

        Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; i++) {
            t->GetEntry(i);
            nTot++;
            int sp = (int)lSp->GetValue();
            if (sp != 0 && sp != 1) { nUnknown++; continue; } // unknown (-1)
            double pt = lPt->GetValue();

            if (sp == 0) nE++; else nP++;
            if (pt < ptLo) ptLo = pt;
            if (pt > ptHi) ptHi = pt;

            // momentum window table and pT plot: tagged tracks with usable layer 0
            if (lQ[0][0]->GetValue() > 0) {
                for (int w = 0; w < nWin; w++)
                    if (pt >= edge[w] && pt < edge[w + 1]) {
                        if (sp == 0) winE[w]++; else winP[w]++;
                    }
                if (sp == 0) { hPtE->Fill(pt); ptSumE += pt; }
                else         { hPtP->Fill(pt); ptSumP += pt; }
            }

            // charge histograms, with the optional pT window
            if (pt < ptMin || pt > ptMax) continue;
            for (int s = 0; s < nSlice; s++) {
                for (int l = 0; l < nLayer; l++) {
                    double q = lQ[s][l]->GetValue();
                    if (q <= 0) continue;
                    if (sp == 0) hE[s][l]->Fill(q);
                    else         hP[s][l]->Fill(q);
                }
            }
        }
    }

    // ---- printout ----
    printf("\n=== Species ===\n"
           "Tracks total : %ld\n"
           "Unknown (-1) : %ld  (%.1f%%)\n"
           "Electrons    : %ld\n"
           "Pions        : %ld\n"
           "fPt range of tagged tracks: %.2f to %.2f GeV/c\n",
           nTot, nUnknown, 100.0 * nUnknown / std::max(nTot, 1L), nE, nP, ptLo, ptHi);

    printf("\n=== Tagged tracks with usable layer 0, by pT ===\n"
           "%-18s %10s %10s\n", "pT window", "electrons", "pions");
    for (int w = 0; w < nWin; w++) {
        if (w < nWin - 1)
            printf("%4.1f - %4.1f GeV/c %10ld %10ld\n", edge[w], edge[w + 1], winE[w], winP[w]);
        else
            printf("%4.1f -  inf GeV/c %10ld %10ld\n", edge[w], winE[w], winP[w]);
    }
    printf("mean pT: electrons %.2f, pions %.2f GeV/c\n",
           hPtE->GetEntries() > 0 ? ptSumE / hPtE->GetEntries() : 0.0,
           hPtP->GetEntries() > 0 ? ptSumP / hPtP->GetEntries() : 0.0);

    printf("\n=== Charge histograms (pT window %.2f to %.2f GeV/c) ===\n", ptMin, ptMax);

    // ---- canvas 1: 3x6 charge grid ----
    TCanvas* c = new TCanvas("c", "TRD charge separation", 2400, 1200);
    c->Divide(nLayer, nSlice);

    for (int s = 0; s < nSlice; s++) {
        for (int l = 0; l < nLayer; l++) {
            c->cd(s * nLayer + l + 1);
            TH1D* e = hE[s][l];
            TH1D* p = hP[s][l];
            printf("Q%d%d: electrons %.0f, pions %.0f\n", s, l, e->GetEntries(), p->GetEntries());
            if (e->Integral() > 0) e->Scale(1.0 / e->Integral());
            if (p->Integral() > 0) p->Scale(1.0 / p->Integral());
            e->SetLineColor(kRed + 1);  e->SetLineWidth(2);
            p->SetLineColor(kBlue + 1); p->SetLineWidth(2);
            double m = std::max(e->GetMaximum(), p->GetMaximum());
            e->SetMaximum(1.2 * m);
            e->Draw("hist");
            p->Draw("hist same");

            TLegend* leg = new TLegend(0.60, 0.72, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.05);
            leg->AddEntry(e, "Electron", "l");
            leg->AddEntry(p, "Pion", "l");
            leg->Draw();
        }
    }
    c->SaveAs("plots/trd_separation03_all_q.png");

    // ---- canvas 2: pT distributions ----
    TCanvas* c2 = new TCanvas("c2", "pT of tagged tracks", 900, 600);
    TH1D* pe = (TH1D*)hPtE->Clone("pe");
    TH1D* pp = (TH1D*)hPtP->Clone("pp");
    if (pe->Integral() > 0) pe->Scale(1.0 / pe->Integral());
    if (pp->Integral() > 0) pp->Scale(1.0 / pp->Integral());
    pe->SetLineColor(kRed + 1);  pe->SetLineWidth(2);
    pp->SetLineColor(kBlue + 1); pp->SetLineWidth(2);
    pe->SetMaximum(1.2 * std::max(pe->GetMaximum(), pp->GetMaximum()));
    pe->Draw("hist");
    pp->Draw("hist same");

    TLegend* leg2 = new TLegend(0.62, 0.75, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(pe, Form("Electron (%.0f)", hPtE->GetEntries()), "l");
    leg2->AddEntry(pp, Form("Pion (%.0f)", hPtP->GetEntries()), "l");
    leg2->Draw();
    c2->SaveAs("plots/trd_separation03_pt.png");
}