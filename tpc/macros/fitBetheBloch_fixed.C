/*
 * File: fitBetheBloch_fixed.C
 * Author: Madhva Fakare
 *
 * Reads TPC dE/dx data from the TPC+TOF skim (O2tpctofskimwde) and fits
 * the Bethe-Bloch formula to extract detector response parameters.
 * This skim contains all tracks regardless of species -- no V0 topology
 * or PID tagging is used here.
 *
 * The pion band dominates the spectrum because the skim is populated
 * mainly by tracks from K0s → π+π− decays, but the histogram is filled
 * with ALL tracks. The pion MPV is extracted by finding the dominant
 * peak in each momentum slice.
 *
 * Theoretical curves for other species are derived by swapping the mass
 * in the same fitted BB parameters. The four species and their masses
 * (from tpcSkimsTableCreator.cxx / PhysicsConstants.h):
 *
 *   Species  | Mass (GeV/c²) | Constant name
 *   ---------+---------------+----------------
 *   e±       | 0.000511      | MassElectron
 *   π±       | 0.13957       | MassPiPlus
 *   K±       | 0.49368       | MassKPlus
 *   p / p̄   | 0.93827       | MassProton
 *
 * === Bethe-Bloch formula (ALICE/ALEPH convention) ===
 *
 *   f(βγ) = P1/β^P4 · (P2 − β^P4 − ln(P3 + (βγ)^(−P5)))
 *
 * The O2Physics implementation is in:
 *   O2Physics/Common/Core/PID/TPCPIDResponse.h
 * which uses BetheBlochAleph from MathUtils with default parameters:
 *
 *   P1 = 0.03209809958934784    (overall amplitude scale)
 *   P2 = 19.9768009185791       (plateau level)
 *   P3 = 2.5266601063857674e-16 (log argument offset, effectively 0)
 *   P4 = 2.7212300300598145     (1/β² exponent, theoretical value = 2.0)
 *   P5 = 6.080920219421387      (relativistic rise exponent)
 *
 * These are the CCDB calibrated values used internally by O2Physics to
 * compute nσ. Note: the O2 implementation multiplies by mMIP = 50.0
 * (arb. units) which sets the absolute signal scale. Our fitted
 * parameters encode the MIP scale directly in P1 instead.
 *
 * The O2 default P4 = 2.72 confirms the physical range for P4 is
 * approximately 2.0-3.0 for the ALICE TPC Ne-CO2-N2 gas mixture.
 *
 * === Fit strategy ===
 *   Stage 1: Fix P4 = 2.0 (theoretical 1/β² exponent) and fit P1, P2,
 *            P3, P5 only. This reduces the 5D parameter space to 4D and
 *            prevents Minuit from wandering into unphysical minima.
 *            Fit range 0.3-4.0 GeV/c: below 0.3 GeV/c the kaon/proton
 *            bands contaminate the pion MPV; above 4.0 GeV/c statistics
 *            are sparse.
 *   Stage 2: Release P4 with limits [0.5, 3.0]. The O2Physics default
 *            P4 = 2.72 confirms this range is physically appropriate for
 *            the ALICE TPC Ne-CO2-N2 gas mixture. P4 converges at ~3.0
 *            which is consistent with the known relativistic correction.
 *
 * === Note on species curves ===
 *   K±, p/p̄, e± curves are NOT separately fitted. The same P1-P5 from
 *   the pion fit are used with only the mass swapped. This is physically
 *   correct: BB parameters describe the detector gas response and are
 *   species-independent. Only the mass sets βγ = p/m at a given momentum.
 *
 * === Note on fPidIndex ===
 *   The TPC+TOF skim (O2tpctofskimwde) does carry an fPidIndex branch,
 *   confirmed via TBrowser inspection. The fPidIndex values follow the
 *   o2::track::PID enum defined in PID.h line 88:
 *
 *   fPidIndex | Species  | Source
 *   ----------+----------+---------------------------
 *       0     | electron | o2::track::PID::Electron
 *       2     | pion     | o2::track::PID::Pion
 *       4     | proton   | o2::track::PID::Proton
 *
 *   PID.h location:
 *     /home/madhvafakare/alice/O2/DataFormats/Reconstruction/include/
 *     ReconstructionDataFormats/PID.h
 *
 *   However, fPidIndex is NOT used in this macro -- all tracks are
 *   filled into one histogram regardless of species. The TBrowser
 *   inspection showed peaks at index 2 (pion) and 4 (proton) in the
 *   TOF skim, confirming pions and protons dominate as expected.
 *   For species-separated histograms using fPidIndex, see learn_bb2.C
 *   and learn_bb3.C which use O2tpcskimv0wde.
 *
 * === Note on species-separated analysis ===
 *   For histograms with individual species isolated via V0 decay topology,
 *   see learn_bb3.C which uses O2tpcskimv0wde with fPidIndex tagging.
 *
 * Run with:
 *   root -l -q macros/fitBetheBloch_fixed.C
 */

Double_t BetheBloch(Double_t* x, Double_t* par) {
    Double_t p    = x[0];
    Double_t mass = par[5];
    Double_t bg   = p / mass;
    Double_t beta = bg / TMath::Sqrt(1.0 + bg * bg);
    if (beta <= 0 || beta >= 1) return 0;
    Double_t P1 = par[0];
    Double_t P2 = par[1];
    Double_t P3 = par[2];
    Double_t P4 = par[3];
    Double_t P5 = par[4];
    Double_t logArg = P3 + TMath::Power(1.0 / bg, P5);
    if (logArg <= 0) return 0;
    return (P1 / TMath::Power(beta, P4)) *
           (P2 - TMath::Power(beta, P4) -
            TMath::Log(logArg));
}

void fitBetheBloch_fixed() {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gRandom->SetSeed(42);

    /* --- Cache: if the filled histogram already exists on disk, load it
     *     directly and skip the ~7-minute data loop entirely.
     *     Delete plots/dedx_cache.root to force a re-fill. -------------- */
    const char* kCacheFile = "plots/dedx_cache.root";
    TH2F* hDeDx = nullptr;

    TFile* fCache = TFile::Open(kCacheFile);
    if (fCache && !fCache->IsZombie()) {
        hDeDx = (TH2F*)fCache->Get("hDeDx");
        if (hDeDx) {
            hDeDx->SetDirectory(0);
            fCache->Close();
            cout << "Loaded histogram from cache (" << kCacheFile << ")." << endl;
            cout << "Total entries: " << hDeDx->GetEntries() << endl;
        } else {
            fCache->Close();
            hDeDx = nullptr;
        }
    }

    /* Always update title and axis labels -- cached histogram keeps old strings */
    if (hDeDx) {
        hDeDx->SetTitle("TPC d#it{E}/d#it{x} vs Momentum - Pb-Pb 2023, LHC23zzf-zzo (apass5);"
                        "#it{p} (GeV/#it{c});TPC Signal (arb. units)");
    }

    if (!hDeDx) {
        hDeDx = new TH2F("hDeDx",
            "TPC d#it{E}/d#it{x} vs Momentum - Pb-Pb 2023, LHC23zzf-zzo (apass5);"
            "#it{p} (GeV/#it{c});TPC Signal (arb. units)",
            100, 0.1, 5, 100, 0, 300);
        hDeDx->SetDirectory(0);

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
            TFile* f = TFile::Open(files[iFile]);
            if (!f || f->IsZombie()) continue;
            cout << "File " << iFile + 1 << "/" << files.size() << endl;

            TIter next(f->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)next())) {
                TString dirName = key->GetName();
                if (!dirName.BeginsWith("DF_")) continue;
                TTree* tree = (TTree*)f->Get(dirName + "/O2tpctofskimwde");
                if (!tree) continue;
                tree->SetBranchStatus("*",              0);
                tree->SetBranchStatus("fTPCSignal",     1);
                tree->SetBranchStatus("fTPCInnerParam", 1);

                Float_t tpcSignal, tpcMom;
                tree->SetBranchAddress("fTPCSignal",     &tpcSignal);
                tree->SetBranchAddress("fTPCInnerParam", &tpcMom);

                Long64_t nEntries = tree->GetEntries();
                for (Long64_t i = 0; i < nEntries; i++) {
                    tree->GetEntry(i);
                    if (tpcMom > 0.1 && tpcMom < 5.0)
                        hDeDx->Fill(tpcMom, tpcSignal);
                }
            }
            f->Close();
        }
        cout << "Total entries: " << hDeDx->GetEntries() << endl;

        TFile* fOut = TFile::Open(kCacheFile, "RECREATE");
        hDeDx->Write();
        fOut->Close();
        cout << "Cache saved to " << kCacheFile << endl;
    }

    /* --- Pion MPV extraction ---
     * Pions (from K0s) are the dominant species, so the global maximum of
     * each momentum slice reliably tracks the pion band.
     * For each slice:
     *   (a) find the global maximum bin as a seed,
     *   (b) fit a Gaussian in a narrow window around that seed to extract
     *       a more precise MPV than the raw bin centre. */
    int nBins = hDeDx->GetNbinsX();
    vector<Double_t> pVec, dedxVec, errVec;

    for (int i = 1; i <= nBins; i++) {
        TH1D* slice = hDeDx->ProjectionY("slice", i, i);
        if (slice->GetEntries() < 200) { delete slice; continue; }

        Double_t p       = hDeDx->GetXaxis()->GetBinCenter(i);
        int      maxBin  = slice->GetMaximumBin();
        Double_t seed    = slice->GetBinCenter(maxBin);
        Double_t seedAmp = slice->GetBinContent(maxBin);

        TF1* gaus = new TF1("gaus", "gaus", seed - 15, seed + 15);
        gaus->SetParameters(seedAmp, seed, 5.0);
        int fitStatus = slice->Fit(gaus, "RQ0");

        if (fitStatus == 0 && gaus->GetParameter(2) > 0.5 && gaus->GetParameter(2) < 30.0) {
            pVec.push_back(p);
            dedxVec.push_back(gaus->GetParameter(1));
            errVec.push_back(gaus->GetParError(1));
        } else {
            pVec.push_back(p);
            dedxVec.push_back(seed);
            errVec.push_back(slice->GetBinWidth(maxBin));
        }
        delete gaus;
        delete slice;
    }
    cout << "Pion MPV points: " << pVec.size() << endl;

    TGraphErrors* gMPV_pi = new TGraphErrors(
        pVec.size(), pVec.data(), dedxVec.data(), nullptr, errVec.data());
    gMPV_pi->SetMarkerStyle(20);
    gMPV_pi->SetMarkerSize(0.8);
    gMPV_pi->SetMarkerColor(kRed);

    /* --- Two-stage pion Bethe-Bloch fit ---
     *
     * Stage 1: Fix P4 = 2.0, fit P1 P2 P3 P5.
     * Stage 2: Release P4 with limits [0.5, 3.0].
     *
     * Starting seeds: P1=50 (MIP scale), P2=2.5, P3=1e-3, P5=2.0.
     * These are empirical seeds that reliably find the correct basin.
     * The O2Physics default P4=2.72 confirms convergence near 3.0
     * is physically reasonable for the ALICE TPC gas mixture. */
    TF1* fitBB_pi = new TF1("fitBB_pi", BetheBloch, 0.1, 5.0, 6);

    fitBB_pi->SetParameters(50.0, 2.5, 1e-3, 2.0, 2.0, 0.13957);
    fitBB_pi->FixParameter(3, 2.0);
    fitBB_pi->FixParameter(5, 0.13957);
    fitBB_pi->SetLineColor(kRed);
    fitBB_pi->SetLineWidth(2);

    cout << "\n--- Stage 1: fit with P4 fixed at 2.0 ---" << endl;
    gMPV_pi->Fit(fitBB_pi, "R", "", 0.3, 4.0);

    fitBB_pi->ReleaseParameter(3);
    fitBB_pi->SetParLimits(3, 0.5, 3.0);
    cout << "\n--- Stage 2: release P4 (limits [0.5, 3.0]) ---" << endl;
    gMPV_pi->Fit(fitBB_pi, "R", "", 0.3, 4.0);

    /* Check convergence */
    bool hitLimit = false;
    for (int i = 0; i < 5; i++) {
        Double_t lo, hi;
        fitBB_pi->GetParLimits(i, lo, hi);
        if (lo != hi && lo != 0 && hi != 0) {
            Double_t val = fitBB_pi->GetParameter(i);
            if (TMath::Abs(val - lo) < 1e-6 || TMath::Abs(val - hi) < 1e-6) {
                cout << "  WARNING: P" << i+1 << " = " << val
                     << " hit limit [" << lo << ", " << hi << "]" << endl;
                hitLimit = true;
            }
        }
    }
    if (!hitLimit) cout << "  All parameters converged within limits." << endl;

    /* --- Theoretical curves for all other species ---
     * Same P1-P5 from the pion fit, only the mass changes. */
    Double_t P1 = fitBB_pi->GetParameter(0);
    Double_t P2 = fitBB_pi->GetParameter(1);
    Double_t P3 = fitBB_pi->GetParameter(2);
    Double_t P4 = fitBB_pi->GetParameter(3);
    Double_t P5 = fitBB_pi->GetParameter(4);

    /* Kaon -- MassKPlus = 0.49368 GeV/c² */
    TF1* curveBB_ka = new TF1("curveBB_ka", BetheBloch, 0.1, 5.0, 6);
    curveBB_ka->SetParameters(P1, P2, P3, P4, P5, 0.49368);
    curveBB_ka->FixParameter(5, 0.49368);
    curveBB_ka->SetLineColor(kBlue + 1);
    curveBB_ka->SetLineWidth(2);
    curveBB_ka->SetLineStyle(2);

    /* Proton -- MassProton = 0.93827 GeV/c² */
    TF1* curveBB_pr = new TF1("curveBB_pr", BetheBloch, 0.1, 5.0, 6);
    curveBB_pr->SetParameters(P1, P2, P3, P4, P5, 0.93827);
    curveBB_pr->FixParameter(5, 0.93827);
    curveBB_pr->SetLineColor(kGreen + 2);
    curveBB_pr->SetLineWidth(2);
    curveBB_pr->SetLineStyle(3);

    /* Electron -- MassElectron = 0.000511 GeV/c²
     * Ultra-relativistic across full TPC range, sits at Fermi plateau
     * slightly above the pion MIP band. */
    TF1* curveBB_el = new TF1("curveBB_el", BetheBloch, 0.01, 5.0, 6);
    curveBB_el->SetParameters(P1, P2, P3, P4, P5, 0.000511);
    curveBB_el->FixParameter(5, 0.000511);
    curveBB_el->SetLineColor(kViolet + 1);
    curveBB_el->SetLineWidth(2);
    curveBB_el->SetLineStyle(4);

    /* --- Diagnostic: dE/dx at 1 GeV/c ---
     * Electron and pion differ by < 1 arb. unit at 1 GeV/c,
     * confirming TPC alone cannot separate them at this momentum.
     * O2Physics default MIP = 50.0 arb. units (TPCPIDResponse.h). */
    cout << "\n=== dE/dx at p = 1 GeV/c ===" << endl;
    cout << "  electron: " << curveBB_el->Eval(1.0) << endl;
    cout << "  pion:     " << fitBB_pi->Eval(1.0)   << endl;
    cout << "  kaon:     " << curveBB_ka->Eval(1.0)  << endl;
    cout << "  proton:   " << curveBB_pr->Eval(1.0)  << endl;
    cout << "  (electron should be slightly above pion)" << endl;

    /* --- Draw --- */
    TCanvas* c1 = new TCanvas("c1", "Bethe-Bloch Fit - PbPb 2023", 900, 700);
    c1->SetLogx();
    c1->SetLogz();
    c1->SetGrid();
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.13);

    hDeDx->Draw("COLZ");
    fitBB_pi->Draw("same");
    curveBB_ka->Draw("same");
    curveBB_pr->Draw("same");
    curveBB_el->Draw("same");

    TLatex* tex = new TLatex();
    tex->SetNDC(false);
    tex->SetTextSize(0.038);
    tex->SetTextFont(42);

    tex->SetTextColor(kViolet + 1);
    Double_t elVal = curveBB_el->Eval(0.15);
    if (elVal > 0 && elVal < 300) {
        tex->DrawLatex(0.13, elVal + 6, "e^{#pm}");
    } else {
        tex->DrawLatex(0.13, fitBB_pi->Eval(0.15) + 12, "e^{#pm}");
    }
    tex->SetTextColor(kRed);
    tex->DrawLatex(3.0, fitBB_pi->Eval(3.0) + 6, "#pi^{#pm}");
    tex->SetTextColor(kBlue + 1);
    tex->DrawLatex(0.62, curveBB_ka->Eval(0.62) + 10, "K^{#pm}");
    tex->SetTextColor(kGreen + 2);
    tex->DrawLatex(1.2, curveBB_pr->Eval(1.2) + 10, "p/#bar{p}");

 /*   cout << "\n=== Pion BB fit parameters (applied to all species) ===" << endl;
    cout << "  Reference O2Physics defaults (TPCPIDResponse.h):" << endl;
    cout << "  P1=0.0321, P2=19.977, P3=2.53e-16, P4=2.721, P5=6.081" << endl;
    cout << "  Fitted values:" << endl; */
    for (int i = 0; i < 5; i++)
        cout << "  P" << i + 1 << " = " << fitBB_pi->GetParameter(i)
             << " +/- " << fitBB_pi->GetParError(i) << endl;

    c1->Update();
    c1->SaveAs("plots/bethebloch_fit_fixed.png");
    c1->SaveAs("plots/bethebloch_fit_fixed.pdf");
    cout << "Saved plots/bethebloch_fit_fixed.png" << endl;
}