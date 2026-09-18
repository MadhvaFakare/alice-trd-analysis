/*
 * File: learn_bb3.C
 * Author: Madhva Fakare
 *
 * Reads TPC dE/dx data from the V0 skim (O2tpcskimv0wde) across all 9
 * LHC23zz runs (apass5) and produces species-separated 2D histograms
 * with Bethe-Bloch curves overlaid.
 *
 * Unlike fitBetheBloch_fixed.C which uses the TPC+TOF skim with all
 * tracks mixed, this macro uses the V0 skim where each track is tagged
 * by its parent V0 decay topology via fPidIndex.
 *
 * PID index assignments (from o2::track::PID, PID.h line 88):
 *
 *   fPidIndex | Species | Source V0             | Mass (GeV/c²)
 *   ----------+---------+-----------------------+--------------
 *       0     | e±      | γ  → e+e-             | 0.000511
 *       2     | π±      | K0s → π+π-            | 0.13957
 *       4     | p/p̄     | Λ  → p + π-           | 0.93827
 *
 * PID.h location:
 *   /home/madhvafakare/alice/O2/DataFormats/Reconstruction/include/
 *   ReconstructionDataFormats/PID.h
 *
 * === BB curve parameters ===
 *   BB curves are NOT fitted here. Parameters come directly from
 *   fitBetheBloch_fixed.C which fitted the ALEPH BB formula to the
 *   pion MPV points extracted from the TPC+TOF skim (O2tpctofskimwde)
 *   using all 9 LHC23zz runs -- 930 million tracks.
 *
 *   Fitted parameters (ALICE/ALEPH convention):
 *     P1 = 3.84314    (amplitude scale, encodes MIP = 50 arb. units)
 *     P2 = 7.55448    (plateau shape)
 *     P3 = 0.000789   (log argument offset)
 *     P4 = 3.0        (1/β² exponent, theoretical value = 2.0)
 *     P5 = 4.0023     (relativistic rise exponent)
 *
 *   Reference O2Physics defaults (TPCPIDResponse.h) for comparison:
 *     P1=0.0321, P2=19.977, P3=2.53e-16, P4=2.721, P5=6.081
 *   (different convention: O2 uses a separate mMIP=50 scale factor)
 *
 *   The same P1-P5 are used for all species -- only the mass is swapped.
 *   This is physically correct: BB parameters describe the detector gas
 *   response which is species-independent. Only the mass sets βγ=p/m.
 *
 * === Method ===
 *   1. Loop over all 9 run files and all DF_ folders within each.
 *   2. Fill separate 2D histograms per species using fPidIndex.
 *      hAll = all tracks, hEl = electrons, hPi = pions, hPr = protons.
 *   3. Cache the filled histograms to plots/bb3_cache.root so the
 *      ~10 minute data loop is skipped on subsequent runs.
 *      Delete the cache file to force a re-fill.
 *   4. Also extract pion and electron MPV points slice by slice via
 *      Gaussian fits -- these are drawn as markers to show where the
 *      band centres actually sit in the V0 skim data.
 *   5. Overlay the BB curves using the parameters from fitBetheBloch_fixed.C.
 *
 * === Output ===
 *   Left pad:  all-species 2D histogram from V0 skim, no curves (raw data)
 *   Right pad: same histogram + MPV markers + BB curves for e, π, p
 *
 * Bethe-Bloch empirical formula (ALICE/ALEPH convention):
 *   f(βγ) = P1/β^P4 · (P2 - β^P4 - ln(P3 + (βγ)^(-P5)))
 *
 * Run with:
 *   root -l -q macros/learn_bb3.C
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

void learn_bb3() {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gRandom->SetSeed(42);

    const UChar_t kPidElectron = 0;  // o2::track::PID::Electron
    const UChar_t kPidPion     = 2;  // o2::track::PID::Pion
    const UChar_t kPidProton   = 4;  // o2::track::PID::Proton

    // ── CACHE ────────────────────────────────────────────────────
    const char* kCache = "plots/bb3_cache.root";

    TH2F* hAll = nullptr;
    TH2F* hEl  = nullptr;
    TH2F* hPi  = nullptr;
    TH2F* hPr  = nullptr;

    TFile* fCache = TFile::Open(kCache);
    if (fCache && !fCache->IsZombie()) {
        hAll = (TH2F*)fCache->Get("hAll"); if (hAll) hAll->SetDirectory(0);
        hEl  = (TH2F*)fCache->Get("hEl");  if (hEl)  hEl->SetDirectory(0);
        hPi  = (TH2F*)fCache->Get("hPi");  if (hPi)  hPi->SetDirectory(0);
        hPr  = (TH2F*)fCache->Get("hPr");  if (hPr)  hPr->SetDirectory(0);
        fCache->Close();
        if (hAll && hEl && hPi && hPr) {
            cout << "Loaded histograms from cache." << endl;
        } else {
            hAll = hEl = hPi = hPr = nullptr;
        }
    }

    if (!hAll) {
        hAll = new TH2F("hAll", "All species;#it{p} (GeV/#it{c});TPC Signal (arb. units)",
                        150, 0.05, 5.0, 200, 0, 300);
        hEl  = new TH2F("hEl",  "Electrons (V0: #gamma#rightarrowe^{+}e^{-});#it{p} (GeV/#it{c});TPC Signal (arb. units)",
                        150, 0.05, 5.0, 200, 0, 300);
        hPi  = new TH2F("hPi",  "Pions (V0: K^{0}_{s}#rightarrow#pi^{+}#pi^{-});#it{p} (GeV/#it{c});TPC Signal (arb. units)",
                        150, 0.05, 5.0, 200, 0, 300);
        hPr  = new TH2F("hPr",  "Proton (V0: #Lambda#rightarrowp#pi^{-});#it{p} (GeV/#it{c});TPC Signal (arb. units)",
                        150, 0.05, 5.0, 200, 0, 300);

        vector<TString> files = {
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzg_apass5_250514/AO2D_merge_LHC23zzg.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzh_apass5_250514/AO2D_merge_LHC23zzh.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzi_apass5_250514/AO2D_merge_LHC23zzi.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzk_apass5_250514/AO2D_merge_LHC23zzk.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzl_apass5_250514/AO2D_merge_LHC23zzl.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzm_apass5_250514/AO2D_merge_LHC23zzm.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzn_apass5_250514/AO2D_merge_LHC23zzn.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzo_apass5_250514/AO2D_merge_LHC23zzo.root"
        };

        for (int iFile = 0; iFile < (int)files.size(); iFile++) {
            cout << "File " << iFile+1 << "/9..." << endl;
            TFile* f = TFile::Open(files[iFile]);
            if (!f || f->IsZombie()) { cout << "Skipping." << endl; continue; }

            TIter next(f->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)next())) {
                TString name = key->GetName();
                if (!name.BeginsWith("DF_")) continue;
                TTree* tree = (TTree*)f->Get(name + "/O2tpcskimv0wde");
                if (!tree) continue;

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
                    if (pid == kPidPion)     hPi->Fill(mom, sig);
                    if (pid == kPidProton)   hPr->Fill(mom, sig);
                }
            }
            f->Close();
        }

        cout << "Pion entries:     " << hPi->GetEntries() << endl;
        cout << "Electron entries: " << hEl->GetEntries() << endl;
        cout << "Proton entries:   " << hPr->GetEntries() << endl;

        TFile* fOut = TFile::Open(kCache, "RECREATE");
        hAll->Write(); hEl->Write(); hPi->Write(); hPr->Write();
        fOut->Close();
        cout << "Cache saved to " << kCache << endl;
    }

    // ── PION MPV EXTRACTION ───────────────────────────────────────
    // Extract MPV points from the V0-tagged pion histogram.
    // These markers show where the pion band actually sits in the
    // V0 skim -- useful to see how well the BB curves from the
    // TOF skim align with the V0 skim data.
    vector<Double_t> pPi, dedxPi, errPi;
    for (int i = 1; i <= hPi->GetNbinsX(); i++) {
        TH1D* slice = hPi->ProjectionY("slicePi", i, i);
        if (slice->GetEntries() < 200) { delete slice; continue; }
        Double_t p       = hPi->GetXaxis()->GetBinCenter(i);
        int      maxBin  = slice->GetMaximumBin();
        Double_t seed    = slice->GetBinCenter(maxBin);
        Double_t seedAmp = slice->GetBinContent(maxBin);
        TF1* g = new TF1("g", "gaus", seed - 15, seed + 15);
        g->SetParameters(seedAmp, seed, 5.0);
        int status = slice->Fit(g, "RQ0");
        if (status == 0 && g->GetParameter(2) > 0.5 && g->GetParameter(2) < 30) {
            pPi.push_back(p);
            dedxPi.push_back(g->GetParameter(1));
            errPi.push_back(g->GetParError(1));
        } else {
            pPi.push_back(p);
            dedxPi.push_back(seed);
            errPi.push_back(slice->GetBinWidth(maxBin));
        }
        delete g; delete slice;
    }
    cout << "Pion MPV points extracted: " << pPi.size() << endl;

    TGraphErrors* gMPV_pi = new TGraphErrors(pPi.size(), pPi.data(), dedxPi.data(), nullptr, errPi.data());
    gMPV_pi->SetMarkerStyle(20); gMPV_pi->SetMarkerSize(0.6); gMPV_pi->SetMarkerColor(kRed);

    // ── ELECTRON MPV EXTRACTION ───────────────────────────────────
    vector<Double_t> pEl, dedxEl, errEl;
    for (int i = 1; i <= hEl->GetNbinsX(); i++) {
        TH1D* slice = hEl->ProjectionY("sliceEl", i, i);
        if (slice->GetEntries() < 50) { delete slice; continue; }
        Double_t p       = hEl->GetXaxis()->GetBinCenter(i);
        int      maxBin  = slice->GetMaximumBin();
        Double_t seed    = slice->GetBinCenter(maxBin);
        Double_t seedAmp = slice->GetBinContent(maxBin);
        TF1* g = new TF1("g", "gaus", seed - 15, seed + 15);
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

    TGraphErrors* gMPV_el = new TGraphErrors(pEl.size(), pEl.data(), dedxEl.data(), nullptr, errEl.data());
    gMPV_el->SetMarkerStyle(20); gMPV_el->SetMarkerSize(0.6); gMPV_el->SetMarkerColor(kViolet+1);

    // ── BB CURVES -- parameters from fitBetheBloch_fixed.C ───────
    // No fitting is done here. Parameters were fitted on the TPC+TOF
    // skim (O2tpctofskimwde) using 930M tracks from all 9 LHC23zz runs.
    // Same P1-P5 for all species, only mass changes.
    const Double_t kP1 = 3.84314;
    const Double_t kP2 = 7.55448;
    const Double_t kP3 = 0.000789308;
    const Double_t kP4 = 3.0;
    const Double_t kP5 = 4.0023;

    TF1* curvePi = new TF1("curvePi", BetheBloch, 0.05, 5.0, 6);
    curvePi->SetParameters(kP1, kP2, kP3, kP4, kP5, 0.13957);
    curvePi->FixParameter(5, 0.13957);
    curvePi->SetLineColor(kRed); curvePi->SetLineWidth(2);

    TF1* curveEl = new TF1("curveEl", BetheBloch, 0.05, 5.0, 6);
    curveEl->SetParameters(kP1, kP2, kP3, kP4, kP5, 0.000511);
    curveEl->FixParameter(5, 0.000511);
    curveEl->SetLineColor(kViolet+1); curveEl->SetLineWidth(2); curveEl->SetLineStyle(2);

    TF1* curvePr = new TF1("curvePr", BetheBloch, 0.05, 5.0, 6);
    curvePr->SetParameters(kP1, kP2, kP3, kP4, kP5, 0.93827);
    curvePr->FixParameter(5, 0.93827);
    curvePr->SetLineColor(kGreen+2); curvePr->SetLineWidth(2); curvePr->SetLineStyle(3);

    // diagnostic
    cout << "\n=== dE/dx at p = 1 GeV/c (fitBetheBloch_fixed parameters) ===" << endl;
    cout << "  electron: " << curveEl->Eval(1.0) << endl;
    cout << "  pion:     " << curvePi->Eval(1.0) << endl;
    cout << "  proton:   " << curvePr->Eval(1.0) << endl;
    cout << "  (electron should be slightly above pion)" << endl;

    // ── DRAW ─────────────────────────────────────────────────────
    // Left pad:  raw V0 skim data, no curves
    // Right pad: same data + MPV markers + BB curves
    TCanvas* c = new TCanvas("c", "TPC dE/dx - V0 skim + BB curves", 1600, 700);
    c->Divide(2, 1);

    // left pad
    c->cd(1);
    gPad->SetLogx(); gPad->SetLogz(); gPad->SetGrid();
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.15);
    hAll->SetTitle("TPC d#it{E}/d#it{x} - Pb-Pb 2023, V0 skim (all runs);"
                   "#it{p} (GeV/#it{c});TPC Signal (arb. units)");
    hAll->Draw("COLZ");
    TLatex* tex0 = new TLatex();
    tex0->SetTextSize(0.042); tex0->SetTextFont(42); tex0->SetNDC(false);
    tex0->SetTextColor(kBlack);
    tex0->DrawLatex(0.12, 285, "No fits");

    // right pad
    c->cd(2);
    gPad->SetLogx(); gPad->SetLogz(); gPad->SetGrid();
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.15);
    TH2F* hAll2 = (TH2F*)hAll->Clone("hAll2");
    hAll2->SetTitle("TPC d#it{E}/d#it{x} - Pb-Pb 2023, V0 skim + BB curves;"
                    "#it{p} (GeV/#it{c});TPC Signal (arb. units)");
    hAll2->Draw("COLZ");
    gMPV_pi->Draw("P same");
    gMPV_el->Draw("P same");
    curvePi->Draw("same");
    curveEl->Draw("same");
    curvePr->Draw("same");

    TLatex* tex = new TLatex();
    tex->SetTextSize(0.042); tex->SetTextFont(42); tex->SetNDC(false);

    tex->SetTextColor(kRed);
    tex->DrawLatex(2.0, curvePi->Eval(2.0) + 3, "#pi^{#pm}");

    tex->SetTextColor(kViolet+1);
    tex->DrawLatex(1.5, curveEl->Eval(1.5) + 4, "e^{#pm}");

    tex->SetTextColor(kGreen+2);
    Double_t prVal = curvePr->Eval(0.4);
    if (prVal > 0 && prVal < 280)
        tex->DrawLatex(0.4, prVal + 6, "p/#bar{p}");
    else
        tex->DrawLatex(0.8, 130, "p/#bar{p}");

    c->SaveAs("plots/step6_bb3.png");
    c->SaveAs("plots/step6_bb3.pdf");
    cout << "\nDone!" << endl;
}