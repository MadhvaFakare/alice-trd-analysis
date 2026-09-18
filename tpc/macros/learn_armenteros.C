/*
 * File: learn_armenteros.C
 * Author: Madhva Fakare
 *
 * Produces the Armenteros-Podolanski plot from the V0 skim
 * (O2tpcskimv0wde) using a single run file (LHC23zzf).
 *
 * The Armenteros-Podolanski plot shows:
 *   x-axis: fAlphaV0 -- longitudinal momentum asymmetry (-1 to +1)
 *   y-axis: fQtV0    -- transverse momentum of daughters (0 to 0.25)
 *
 * Three distinct populations appear:
 *   - Gamma conversions: small oval near origin (alpha~0, Qt~0)
 *   - K0s:              large symmetric oval (alpha~0, Qt~0.2)
 *   - Lambda/AntiLambda: two asymmetric arms at +/- alpha
 *
 * Left pad:  all V0 tracks -- full picture
 * Right pad: only electrons (fPidIndex=0) -- where your sample sits
 *
 * Run with:
 *   root -l -q macros/learn_armenteros.C
 */

void learn_armenteros() {

    // ── STYLE ─────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gRandom->SetSeed(42);

    // ── PID INDEX ─────────────────────────────────────────────────
    const UChar_t kPidElectron = 0;

    // ── HISTOGRAMS ────────────────────────────────────────────────
    // hArm:    all V0 tracks -- full Armenteros plot
    // hArmEl:  electrons only -- where fPidIndex=0 sits
    TH2F* hArm = new TH2F("hArm",
        "Armenteros-Podolanski - All V0 tracks;"
        "#alpha_{V0} = (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-});"
        "#it{q}_{T} (GeV/#it{c})",
        200, -1.0, 1.0, 200, 0.0, 0.25);

    TH2F* hArmEl = new TH2F("hArmEl",
        "Armenteros-Podolanski - V0 electrons (fPidIndex=0);"
        "#alpha_{V0} = (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-});"
        "#it{q}_{T} (GeV/#it{c})",
        200, -1.0, 1.0, 200, 0.0, 0.25);

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
        tree->SetBranchStatus("*",          0);
        tree->SetBranchStatus("fAlphaV0",   1);
        tree->SetBranchStatus("fQtV0",      1);
        tree->SetBranchStatus("fPidIndex",  1);

        Float_t alpha, qt;
        UChar_t pid;
        tree->SetBranchAddress("fAlphaV0",  &alpha);
        tree->SetBranchAddress("fQtV0",     &qt);
        tree->SetBranchAddress("fPidIndex", &pid);

        for (Long64_t i = 0; i < tree->GetEntries(); i++) {
            tree->GetEntry(i);
            hArm->Fill(alpha, qt);
            if (pid == kPidElectron) hArmEl->Fill(alpha, qt);
        }
    }

    f->Close();
    cout << "DF_ folders processed: " << nFolders << endl;
    cout << "All V0 entries:      " << hArm->GetEntries() << endl;
    cout << "Electron entries:    " << hArmEl->GetEntries() << endl;

    // ── DRAW ──────────────────────────────────────────────────────
    TCanvas* c = new TCanvas("c", "Armenteros-Podolanski - LHC23zzf",
                              1600, 700);
    c->Divide(2, 1);

    // left pad: all V0 tracks
    c->cd(1);
    gPad->SetLogz();
    gPad->SetGrid();
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.15);
    hArm->Draw("COLZ");

    // label the three populations
    TLatex* tex = new TLatex();
    tex->SetTextSize(0.040);
    tex->SetTextFont(42);
    tex->SetNDC(false);

    tex->SetTextColor(kWhite);
    tex->DrawLatex(0.05, 0.20, "K^{0}_{s}");
    tex->DrawLatex(0.55, 0.20, "K^{0}_{s}");
    tex->DrawLatex(0.70, 0.13, "#Lambda");
    tex->DrawLatex(-0.85, 0.13, "#bar{#Lambda}");
    tex->DrawLatex(0.05, 0.04, "#gamma#rightarrowe^{+}e^{-}");

    tex->SetNDC(true);
    tex->SetTextSize(0.035);
    tex->SetTextColor(kGray+2);
    tex->DrawLatex(0.15, 0.88, "Pb-Pb 2023, LHC23zzf apass5");
    tex->DrawLatex(0.15, 0.83, "V0 skim, all species");

    // right pad: electrons only
    c->cd(2);
    gPad->SetLogz();
    gPad->SetGrid();
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.15);
    hArmEl->Draw("COLZ");

    TLatex* tex2 = new TLatex();
    tex2->SetTextSize(0.040);
    tex2->SetTextFont(42);
    tex2->SetNDC(false);
    tex2->SetTextColor(kWhite);
    tex2->DrawLatex(0.05, 0.04, "#gamma#rightarrowe^{+}e^{-}");

    tex2->SetNDC(true);
    tex2->SetTextSize(0.035);
    tex2->SetTextColor(kGray+2);
    tex2->DrawLatex(0.65, 0.88, "Pb-Pb 2023, LHC23zzf apass5");
    tex2->DrawLatex(0.65, 0.83, "V0 electrons only");

    c->SaveAs("plots/step9_armenteros.png");
    c->SaveAs("plots/step9_armenteros.pdf");
    cout << "\nDone! Saved plots/step9_armenteros.png" << endl;
}