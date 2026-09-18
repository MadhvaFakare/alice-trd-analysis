/*
 * File: learn_nsigma.C
 * Author: Madhva Fakare
 *
 * Computes TPC nSigma separation between electrons and pions
 * using the V0-tagged species histograms from bb3_cache.root.
 *
 * Method:
 *   For each momentum slice, fit a Gaussian to the pion dE/dx band
 *   to extract mean (mu_pi) and width (sigma_pi). These define the
 *   TPC resolution at that momentum. Then for each electron track,
 *   compute:
 *
 *     nSigma_e = (dEdx_track - mu_pi) / sigma_pi
 *
 *   A well-separated electron should sit at positive nSigma.
 *   At 1 GeV/c the electron and pion bands nearly overlap, so
 *   nSigma_e ~ 0 -- this is the key physics result motivating TRD.
 *
 * Output:
 *   Left pad:  nSigma_electron vs p for V0 electrons (signal)
 *   Right pad: nSigma_electron vs p for V0 pions (background/reference)
 *   Also: 1D nSigma projection at p ~ 1 GeV/c showing overlap
 *
 * Run with:
 *   root -l -q macros/learn_nsigma.C
 */

void learn_nsigma() {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gRandom->SetSeed(42);

    // ── LOAD CACHE ───────────────────────────────────────────────
    const char* kCache = "plots/bb3_cache.root";
    TFile* fCache = TFile::Open(kCache);
    if (!fCache || fCache->IsZombie()) {
        cout << "ERROR: cache not found at " << kCache << endl;
        cout << "Run learn_bb3.C first to generate the cache." << endl;
        return;
    }

    TH2F* hEl = (TH2F*)fCache->Get("hEl"); if (hEl) hEl->SetDirectory(0);
    TH2F* hPi = (TH2F*)fCache->Get("hPi"); if (hPi) hPi->SetDirectory(0);
    fCache->Close();

    if (!hEl || !hPi) {
        cout << "ERROR: could not load hEl or hPi from cache." << endl;
        return;
    }
    cout << "Loaded hEl and hPi from cache." << endl;

    // ── OUTPUT HISTOGRAMS ─────────────────────────────────────────
    // nSigma range: -5 to +10 (electrons sit above pions so extend upward)
    int   nBinsP    = hPi->GetNbinsX();
    double pLo      = hPi->GetXaxis()->GetXmin();
    double pHi      = hPi->GetXaxis()->GetXmax();

    TH2F* hNsigEl = new TH2F("hNsigEl",
        "n#sigma_{e} (V0 electrons);#it{p} (GeV/#it{c});n#sigma_{e}^{TPC}",
        nBinsP, pLo, pHi, 100, -5.0, 10.0);

    TH2F* hNsigPi = new TH2F("hNsigPi",
        "n#sigma_{e} (V0 pions);#it{p} (GeV/#it{c});n#sigma_{e}^{TPC}",
        nBinsP, pLo, pHi, 100, -5.0, 10.0);

    // ── SLICE-BY-SLICE: FIT PION BAND, FILL NSIGMA HISTS ─────────
    // For each momentum bin:
    //   1. Fit Gaussian to pion slice -> get mu_pi, sigma_pi
    //   2. Loop over electron slice bins, compute nSigma, fill hNsigEl
    //   3. Loop over pion slice bins, compute nSigma, fill hNsigPi

    int nSlicesOk = 0;

    for (int i = 1; i <= nBinsP; i++) {

        double pCenter = hPi->GetXaxis()->GetBinCenter(i);

        // --- fit pion band in this slice ---
        TH1D* slicePi = hPi->ProjectionY(Form("slicePi_%d", i), i, i);
        if (slicePi->GetEntries() < 200) { delete slicePi; continue; }

        int      maxBin  = slicePi->GetMaximumBin();
        double   seed    = slicePi->GetBinCenter(maxBin);
        double   seedAmp = slicePi->GetBinContent(maxBin);

        TF1* gPi = new TF1(Form("gPi_%d", i), "gaus", seed - 20, seed + 20);
        gPi->SetParameters(seedAmp, seed, 6.0);
        int status = slicePi->Fit(gPi, "RQ0");

        double mu_pi    = gPi->GetParameter(1);
        double sigma_pi = gPi->GetParameter(2);

        // require sensible fit
        if (status != 0 || sigma_pi < 0.5 || sigma_pi > 30.0) {
            delete gPi; delete slicePi; continue;
        }

        nSlicesOk++;

        // --- fill nSigma for electrons in this slice ---
        TH1D* sliceEl = hEl->ProjectionY(Form("sliceEl_%d", i), i, i);
        for (int j = 1; j <= sliceEl->GetNbinsX(); j++) {
            double dEdx    = sliceEl->GetBinCenter(j);
            double nEvents = sliceEl->GetBinContent(j);
            if (nEvents < 1) continue;
            double nsig = (dEdx - mu_pi) / sigma_pi;
            // fill nEvents times (weighted fill)
            hNsigEl->Fill(pCenter, nsig, nEvents);
        }

        // --- fill nSigma for pions in this slice (reference) ---
        for (int j = 1; j <= slicePi->GetNbinsX(); j++) {
            double dEdx    = slicePi->GetBinCenter(j);
            double nEvents = slicePi->GetBinContent(j);
            if (nEvents < 1) continue;
            double nsig = (dEdx - mu_pi) / sigma_pi;
            hNsigPi->Fill(pCenter, nsig, nEvents);
        }

        delete gPi; delete slicePi; delete sliceEl;
    }

    cout << "Momentum slices with good pion fit: " << nSlicesOk << endl;

    // ── 1D PROJECTION AT ~1 GEV/C ────────────────────────────────
    // Find the bin closest to p = 1.0 GeV/c
    int binAt1 = hNsigEl->GetXaxis()->FindBin(1.0);
    // use a small window: 0.9 - 1.1 GeV/c
    int binLo  = hNsigEl->GetXaxis()->FindBin(0.9);
    int binHi  = hNsigEl->GetXaxis()->FindBin(1.1);

    TH1D* proj1El = hNsigEl->ProjectionY("proj1El", binLo, binHi);
    TH1D* proj1Pi = hNsigPi->ProjectionY("proj1Pi", binLo, binHi);

    proj1El->SetLineColor(kViolet+1); proj1El->SetLineWidth(2);
    proj1Pi->SetLineColor(kRed);      proj1Pi->SetLineWidth(2);

    // normalise to unit area for shape comparison
    if (proj1El->Integral() > 0) proj1El->Scale(1.0 / proj1El->Integral());
    if (proj1Pi->Integral() > 0) proj1Pi->Scale(1.0 / proj1Pi->Integral());

    // ── DRAW ─────────────────────────────────────────────────────
    TCanvas* c = new TCanvas("c", "TPC nSigma - electron vs pion", 2100, 700);
    c->Divide(3, 1);

    // pad 1: nSigma for V0 electrons
    c->cd(1);
    gPad->SetLogx(); gPad->SetLogz(); gPad->SetGrid();
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.15);
    hNsigEl->SetTitle(
        "TPC n#sigma_{e} - V0 electrons;#it{p} (GeV/#it{c});n#sigma_{e}^{TPC}");
    hNsigEl->Draw("COLZ");
    // draw reference lines at 0 and +2
    TLine* lZero = new TLine(pLo, 0, pHi, 0);
    lZero->SetLineColor(kWhite); lZero->SetLineStyle(2); lZero->SetLineWidth(1);
    lZero->Draw();

    // pad 2: nSigma for V0 pions (shows where the background sits)
    c->cd(2);
    gPad->SetLogx(); gPad->SetLogz(); gPad->SetGrid();
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.15);
    hNsigPi->SetTitle(
        "TPC n#sigma_{e} - V0 pions;#it{p} (GeV/#it{c});n#sigma_{e}^{TPC}");
    hNsigPi->Draw("COLZ");
    TLine* lZero2 = new TLine(pLo, 0, pHi, 0);
    lZero2->SetLineColor(kWhite); lZero2->SetLineStyle(2); lZero2->SetLineWidth(1);
    lZero2->Draw();

    // pad 3: 1D projection at 1 GeV/c - the key plot
    c->cd(3);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.05);

    double yMax = max(proj1El->GetMaximum(), proj1Pi->GetMaximum()) * 1.3;
    proj1El->SetMaximum(yMax);
    proj1El->SetTitle(
        "n#sigma_{e}^{TPC} at #it{p} = 0.9-1.1 GeV/#it{c};"
        "n#sigma_{e}^{TPC};Normalised counts");
    proj1El->Draw("hist");
    proj1Pi->Draw("hist same");

    TLegend* leg = new TLegend(0.55, 0.70, 0.92, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
    leg->AddEntry(proj1El, "V0 electrons (signal)", "l");
    leg->AddEntry(proj1Pi, "V0 pions (background)", "l");
    leg->Draw();

    // vertical line at nSigma = 0
    TLine* lRef = new TLine(0, 0, 0, yMax * 0.95);
    lRef->SetLineColor(kGray+1); lRef->SetLineStyle(2); lRef->SetLineWidth(1);
    lRef->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextSize(0.038); tex->SetTextFont(42); tex->SetNDC(true);
    tex->SetTextColor(kGray+2);
    tex->DrawLatex(0.15, 0.85, "Pb-Pb 2023, LHC23zz apass5");
    tex->DrawLatex(0.15, 0.80, "9 runs, V0 skim");

    c->SaveAs("plots/step7_nsigma.png");
    c->SaveAs("plots/step7_nsigma.pdf");
    cout << "\nDone! Saved plots/step7_nsigma.png" << endl;
}