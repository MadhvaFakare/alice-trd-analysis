/*
 * File: learn_electron.C
 * Author: Madhva Fakare
 *
 * Plots all species in the background from the V0 skim
 * (O2tpcskimv0wde) with the electron BB curve and MPV
 * markers highlighted on top.
 *
 * Goal: show where the electron band sits among all other
 * particles and identify momentum regions where electrons
 * are clearly separable.
 *
 * PID index assignments (from o2::track::PID, PID.h):
 *   fPidIndex = 0 -> electron (gamma -> e+e-)
 *
 * BB parameters from fitBetheBloch_fixed.C (fitted on TOF skim,
 * 930M tracks, all 9 LHC23zz apass5 runs):
 *   P1=3.84314, P2=7.55448, P3=0.000789308, P4=3.0, P5=4.0023
 *
 * Run with:
 *   root -l -q macros/learn_electron.C
 */

Double_t BetheBloch(Double_t* x, Double_t* par) {
    Double_t p    = x[0];
    Double_t mass = par[5];
    Double_t bg   = p / mass;
    Double_t beta = bg / TMath::Sqrt(1.0 + bg * bg);
    if (beta <= 0 || beta >= 1) return 0;
    Double_t betaP4  = TMath::Power(beta, par[3]);
    Double_t invbgP5 = TMath::Power(1.0 / bg, par[4]);
    Double_t logArg  = par[2] + invbgP5;
    if (logArg <= 0) return 0;
    Double_t val = (par[0] / betaP4) * (par[1] - betaP4 - TMath::Log(logArg));
    if (val <= 0) return 0;
    return val;
}

void learn_electron() {

    // ── STYLE ─────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gRandom->SetSeed(42);

    // ── PID INDEX ─────────────────────────────────────────────────
    const UChar_t kPidElectron = 0;

    // ── HISTOGRAMS ────────────────────────────────────────────────
    // hAll: all tracks -- background context
    // hEl:  electrons only -- for MPV extraction
    TH2F* hAll = new TH2F("hAll",
        "All species;"
        "#it{p} (GeV/#it{c});TPC Signal (arb. units)",
        150, 0.05, 5.0, 200, 0, 300);

    TH2F* hEl = new TH2F("hEl",
        "Electrons (V0: #gamma#rightarrowe^{+}e^{-});"
        "#it{p} (GeV/#it{c});TPC Signal (arb. units)",
        150, 0.05, 5.0, 200, 0, 300);

    // ── OPEN FILE ─────────────────────────────────────────────────
    TFile* f = TFile::Open(
        "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/"
        "LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    if (!f || f->IsZombie()) {
        cout << "ERROR: could not open file." << endl;
        return;
    }
    cout << "File opened successfully." << endl;

    // ── LOOP OVER DF_ FOLDERS ─────────────────────────────────────
    TIter next(f->GetListOfKeys());
    TKey* key;
    int nFolders = 0;

    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;

        TTree* tree = (TTree*)f->Get(name + "/O2tpcskimv0wde");
        if (!tree) continue;
        nFolders++;

        // only read the branches we need
        tree->SetBranchStatus("*",              0);
        tree->SetBranchStatus("fTPCSignal",     1);
        tree->SetBranchStatus("fTPCInnerParam", 1);
        tree->SetBranchStatus("fPidIndex",      1);

        Float_t sig, mom;
        UChar_t pid;
        tree->SetBranchAddress("fTPCSignal",     &sig);
        tree->SetBranchAddress("fTPCInnerParam", &mom);
        tree->SetBranchAddress("fPidIndex",      &pid);

        for (Long64_t i = 0; i < tree->GetEntries(); i++) {
            tree->GetEntry(i);
            if (mom < 0.05 || mom > 5.0) continue;
            hAll->Fill(mom, sig);
            if (pid == kPidElectron) hEl->Fill(mom, sig);
        }
    }

    f->Close();
    cout << "DF_ folders processed: " << nFolders << endl;
    cout << "All entries:      " << hAll->GetEntries() << endl;
    cout << "Electron entries: " << hEl->GetEntries() << endl;

    // ── ELECTRON MPV EXTRACTION ───────────────────────────────────
    // Fit Gaussian to each momentum slice of the electron histogram
    // to extract the peak position (MPV) and width (resolution)
    vector<Double_t> pEl, dedxEl, errEl;

    for (int i = 1; i <= hEl->GetNbinsX(); i++) {
        TH1D* slice = hEl->ProjectionY(Form("sliceEl_%d", i), i, i);
        if (slice->GetEntries() < 50) { delete slice; continue; }

        Double_t p       = hEl->GetXaxis()->GetBinCenter(i);
        int      maxBin  = slice->GetMaximumBin();
        Double_t seed    = slice->GetBinCenter(maxBin);
        Double_t seedAmp = slice->GetBinContent(maxBin);

        TF1* g = new TF1(Form("gEl_%d", i), "gaus", seed - 15, seed + 15);
        g->SetParameters(seedAmp, seed, 5.0);
        int status = slice->Fit(g, "RQ0");

        if (status == 0 && g->GetParameter(2) > 0.5 && g->GetParameter(2) < 30) {
            pEl.push_back(p);
            dedxEl.push_back(g->GetParameter(1));
            errEl.push_back(g->GetParError(1));
        } else {
            pEl.push_back(p);
            dedxEl.push_back(seed);
            errEl.push_back(slice->GetBinWidth(maxBin));
        }
        delete g; delete slice;
    }
    cout << "Electron MPV points extracted: " << pEl.size() << endl;

    TGraphErrors* gMPV_el = new TGraphErrors(
        pEl.size(), pEl.data(), dedxEl.data(), nullptr, errEl.data());
    gMPV_el->SetMarkerStyle(20);
    gMPV_el->SetMarkerSize(0.6);
    gMPV_el->SetMarkerColor(kViolet+1);

    // ── BB CURVE ──────────────────────────────────────────────────
    // Parameters from fitBetheBloch_fixed.C
    const Double_t kP1 = 3.84314;
    const Double_t kP2 = 7.55448;
    const Double_t kP3 = 0.000789308;
    const Double_t kP4 = 3.0;
    const Double_t kP5 = 4.0023;

    TF1* curveEl = new TF1("curveEl", BetheBloch, 0.05, 5.0, 6);
    curveEl->SetParameters(kP1, kP2, kP3, kP4, kP5, 0.000511);
    curveEl->FixParameter(5, 0.000511);
    curveEl->SetLineColor(kViolet+1);
    curveEl->SetLineWidth(2);
    curveEl->SetLineStyle(2);

    // diagnostic
    cout << "\n=== dE/dx at p = 1 GeV/c ===" << endl;
    cout << "  electron: " << curveEl->Eval(1.0) << endl;

    // ── DRAW ──────────────────────────────────────────────────────
    TCanvas* c = new TCanvas("c", "Electron band - LHC23zzf", 900, 700);
    c->SetLogx();
    c->SetLogz();
    c->SetGrid();
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.15);

    hAll->SetTitle(
        "TPC d#it{E}/d#it{x} - Pb-Pb 2023 LHC23zzf, electron highlighted;"
        "#it{p} (GeV/#it{c});TPC Signal (arb. units)");
    hAll->Draw("COLZ");

    // electron MPV markers and BB curve on top
    gMPV_el->Draw("P same");
    curveEl->Draw("same");

    // label
    TLatex* tex = new TLatex();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetNDC(false);
    tex->SetTextColor(kViolet+1);
    tex->DrawLatex(1.5, curveEl->Eval(1.5) + 5, "e^{#pm}");

    tex->SetNDC(true);
    tex->SetTextSize(0.035);
    tex->SetTextColor(kGray+2);
    tex->DrawLatex(0.15, 0.88, "Pb-Pb 2023, LHC23zzf apass5");
    tex->DrawLatex(0.15, 0.83, "V0 skim, single run");

    c->SaveAs("plots/step8_electron.png");
    c->SaveAs("plots/step8_electron.pdf");
    cout << "\nDone! Saved plots/step8_electron.png" << endl;
}