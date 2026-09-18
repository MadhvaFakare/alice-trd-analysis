// learn_gammatag.C
//
// Purpose: fMass turned out to be a per-track PID-hypothesis mass, not V0
// invariant mass, so it cannot tag photon conversions. This macro moves to
// tools that are actually built for the job:
//   1) fPidIndex   - does an electron hypothesis even exist in this tree?
//   2) kinematic window (alpha, qT) - define a gamma-candidate region with
//      no mass assumption at all
//   3) fGammaPsiPair - the purpose-built ALICE conversion-tagging variable,
//      check whether it actually discriminates inside that region
//   4) fRadiusV0   - cross-check candidates against detector material radii
//
// This is exploratory: ranges/bins below are starting guesses, widen or
// narrow them based on what you actually see.
//
// Run with: root -l learn_gammatag.C

void learn_gammatag() {

    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    TTree *tree = (TTree*)f->Get("DF_2336518081137728/O2tpcskimv0wde");

    if (!tree) {
        printf("ERROR: could not find tree, check the path\n");
        return;
    }

    Long64_t nTotal = tree->GetEntries();
    printf("Total entries: %lld\n\n", nTotal);

    // ---------------------------------------------------------------
    // STEP 1: what PID hypotheses actually exist in this tree?
    // ---------------------------------------------------------------
    TCanvas *cPid = new TCanvas("cPid", "fPidIndex", 800, 600);
    tree->Draw("fPidIndex>>hPid(20,-0.5,19.5)");
    TH1F *hPid = (TH1F*)gDirectory->Get("hPid");
    hPid->SetTitle("fPidIndex;PID index;Counts");
    cPid->SaveAs("pidindex_check.png");

    printf("--- fPidIndex populated bins ---\n");
    for (int i = 1; i <= hPid->GetNbinsX(); i++) {
        double content = hPid->GetBinContent(i);
        if (content > 0) {
            printf("  index %.0f : %.0f entries\n", hPid->GetBinCenter(i), content);
        }
    }
    printf("Check this list against your PID index convention. If no populated\n");
    printf("index corresponds to electron, fMass/fPidIndex was never going to\n");
    printf("show a gamma-related peak here, regardless of cuts.\n\n");

    // ---------------------------------------------------------------
    // STEP 2: kinematic gamma-candidate region, no mass involved.
    // Uses the actual published ALICE PCM selection: an elliptical
    // boundary in (alpha, qT), not a plain box. A narrow |alpha| box
    // (an earlier version of this macro used |alpha|<0.2) cuts away
    // real conversions at higher |alpha| -- equal-mass two-body-like
    // kinematics does NOT mean alpha stays near zero, K0s in your own
    // plot is proof of that.
    // ---------------------------------------------------------------
    double qTmax = 0.05;
    double alphaMax = 0.95;
    TCut gammaKinematic = Form("TMath::Abs(fAlphaV0)<%f && fQtV0 < %f*TMath::Sqrt(1 - (fAlphaV0*fAlphaV0)/(%f*%f))",
                                alphaMax, qTmax, alphaMax, alphaMax);
    Long64_t nGammaKinematic = tree->GetEntries(gammaKinematic);
    printf("--- Kinematic gamma-candidate region (ALICE PCM ellipse: qTmax=%.2f, alphaMax=%.2f) ---\n", qTmax, alphaMax);
    printf("Entries: %lld  (fraction of total: %.6f)\n\n", nGammaKinematic, (double)nGammaKinematic/nTotal);

    // ---------------------------------------------------------------
    // STEP 3: does fGammaPsiPair actually discriminate in this region?
    // ---------------------------------------------------------------
    TCanvas *cPsi = new TCanvas("cPsi", "fGammaPsiPair", 800, 600);
    tree->Draw("fGammaPsiPair>>hPsiAll(200,-2,2)");
    TH1F *hPsiAll = (TH1F*)gDirectory->Get("hPsiAll");
    hPsiAll->SetLineColor(kBlack);
    hPsiAll->SetTitle("fGammaPsiPair;Psi pair;Counts");

    tree->SetLineColor(kMagenta);
    tree->Draw("fGammaPsiPair>>hPsiKinGamma(200,-2,2)", gammaKinematic, "same");
    TH1F *hPsiKinGamma = (TH1F*)gDirectory->Get("hPsiKinGamma");
    hPsiKinGamma->SetLineColor(kMagenta);

    cPsi->SaveAs("psipair_check.png");

    printf("--- fGammaPsiPair check ---\n");
    printf("Black = full sample. Magenta = kinematic gamma-candidate subset only.\n");
    printf("If magenta is visibly narrower/more peaked than black, PsiPair is\n");
    printf("doing real discriminating work and is usable as a cut.\n");
    printf("If the range above looks clipped (population piling at -2 or +2),\n");
    printf("widen the histogram range and rerun before concluding anything.\n\n");

    // ---------------------------------------------------------------
    // STEP 4: radius cross-check, independent of PID and mass entirely
    // ---------------------------------------------------------------
    TCanvas *cRad = new TCanvas("cRad", "fRadiusV0 comparison", 800, 600);
    tree->SetLineColor(kBlack);
    tree->Draw("fRadiusV0>>hRadAll(200,0,200)");
    TH1F *hRadAll = (TH1F*)gDirectory->Get("hRadAll");
    hRadAll->SetTitle("fRadiusV0;Radius (cm);Counts");

    tree->SetLineColor(kMagenta);
    tree->Draw("fRadiusV0>>hRadGamma(200,0,200)", gammaKinematic, "same");
    TH1F *hRadGamma = (TH1F*)gDirectory->Get("hRadGamma");
    hRadGamma->SetLineColor(kMagenta);

    cRad->SaveAs("radius_check.png");

    printf("--- fRadiusV0 check ---\n");
    printf("Black = full sample. Magenta = kinematic gamma-candidate subset.\n");
    printf("Real conversions cluster at detector-material radii (beampipe,\n");
    printf("ITS layers). Structure in magenta that black does not share\n");
    printf("supports a genuine conversion population, independent of PID/mass.\n\n");

    // ---------------------------------------------------------------
    // STEP 5: take the corrected (alpha, qT) selection back to the
    // actual Armenteros plot this was all in service of
    // ---------------------------------------------------------------
    TCanvas *cArmRefined = new TCanvas("cArmRefined", "Armenteros with ALICE PCM gamma selection", 1000, 800);
    tree->Draw("fQtV0:fAlphaV0>>hAllRefined(200,-1,1,200,0,0.3)", "", "colz");
    TH2F *hAllRefined = (TH2F*)gDirectory->Get("hAllRefined");
    hAllRefined->SetTitle("Armenteros-Podolanski with ALICE PCM gamma selection;#alpha;q_{T} (GeV/c)");

    tree->SetMarkerColor(kMagenta);
    tree->SetMarkerStyle(20);
    tree->SetMarkerSize(0.4);
    tree->Draw("fQtV0:fAlphaV0", gammaKinematic, "same");

    cArmRefined->SaveAs("armenteros_refined_gammatag.png");

    // optional cross-check against the material-radius clustering found
    // earlier -- extra, not part of the published ALICE selection itself
    TCut gammaRadius = "fRadiusV0>15 && fRadiusV0<35";
    Long64_t nBoth = tree->GetEntries(gammaKinematic && gammaRadius);
    printf("--- Interpretation ---\n");
    printf("This should now show magenta points spread across a WIDE alpha\n");
    printf("range (out toward +-0.9) while staying at low qT, forming a flat\n");
    printf("band under the K0s arc, not a narrow blob at alpha=0.\n");
    printf("Of the %lld selected candidates, %lld (%.1f%%) also fall in the\n",
           nGammaKinematic, nBoth, 100.0*nBoth/nGammaKinematic);
    printf("15-35 cm radius window seen earlier, a further, non-standard\n");
    printf("purity cross-check, not required by the published selection.\n");
    printf("Also worth a visual check: at higher |alpha| (~0.6-0.8), does the\n");
    printf("magenta selection creep into the very bottom tips of the Lambda/\n");
    printf("AntiLambda horns? The ellipse's allowed qT there is close to the\n");
    printf("horns' minimum qT, so some edge overlap is possible.\n");

    // ---------------------------------------------------------------
    // STEP 6: test the cheapest fix for horn-tip leakage: does requiring
    // fPidIndex==0 (electron-like dE/dx) on top of the ellipse remove it?
    // A genuine Lambda proton daughter should be nowhere near electron-like
    // dE/dx, unlike the earlier pion/electron overlap problem.
    // ---------------------------------------------------------------
    TCut gammaPlusPid = gammaKinematic && "fPidIndex==0";
    Long64_t nGammaPlusPid = tree->GetEntries(gammaPlusPid);
    printf("\n--- Ellipse + fPidIndex==0 ---\n");
    printf("Entries: %lld  (was %lld before this extra cut)\n\n", nGammaPlusPid, nGammaKinematic);

    TCanvas *cArmPid = new TCanvas("cArmPid", "Armenteros with ellipse + PID cut", 1000, 800);
    tree->Draw("fQtV0:fAlphaV0>>hAllPid(200,-1,1,200,0,0.3)", "", "colz");
    TH2F *hAllPid = (TH2F*)gDirectory->Get("hAllPid");
    hAllPid->SetTitle("Armenteros-Podolanski, ellipse + fPidIndex==0;#alpha;q_{T} (GeV/c)");

    tree->SetMarkerColor(kMagenta);
    tree->SetMarkerStyle(20);
    tree->SetMarkerSize(0.4);
    tree->Draw("fQtV0:fAlphaV0", gammaPlusPid, "same");
    cArmPid->SaveAs("armenteros_ellipse_plus_pid.png");

    printf("Check: do the two blocks sitting on the Lambda/AntiLambda horn\n");
    printf("tips disappear now? If yes, fPidIndex==0 is doing real work.\n\n");

    // ---------------------------------------------------------------
    // STEP 7: understand what fNSigTPC actually represents before relying
    // on it further -- is it a fixed hypothesis, or "sigma relative to
    // whatever fPidIndex already picked" (which would trivially cluster
    // near 0 for every group and tell us nothing new)?
    // ---------------------------------------------------------------
    TCanvas *cNSig = new TCanvas("cNSig", "fNSigTPC by fPidIndex group", 800, 600);
    tree->SetLineColor(kRed);
    tree->Draw("fNSigTPC>>hNSig0(200,-10,10)", "fPidIndex==0");
    TH1F *hNSig0 = (TH1F*)gDirectory->Get("hNSig0");
    hNSig0->SetLineColor(kRed);
    hNSig0->SetTitle("fNSigTPC by fPidIndex group (red=0, blue=2, green=4);n#sigma_{TPC};Counts");

    tree->SetLineColor(kBlue);
    tree->Draw("fNSigTPC>>hNSig2(200,-10,10)", "fPidIndex==2", "same");
    TH1F *hNSig2 = (TH1F*)gDirectory->Get("hNSig2");
    hNSig2->SetLineColor(kBlue);

    tree->SetLineColor(kGreen+2);
    tree->Draw("fNSigTPC>>hNSig4(200,-10,10)", "fPidIndex==4", "same");
    TH1F *hNSig4 = (TH1F*)gDirectory->Get("hNSig4");
    hNSig4->SetLineColor(kGreen+2);

    cNSig->SaveAs("nsigtpc_by_pidindex.png");

    printf("--- fNSigTPC by fPidIndex group ---\n");
    printf("If all three curves center near 0, fNSigTPC is relative to\n");
    printf("whatever fPidIndex already chose (self-referential, not directly\n");
    printf("useful as an independent cut). If they sit at DIFFERENT centers,\n");
    printf("fNSigTPC is a fixed hypothesis (likely electron), usable directly,\n");
    printf("e.g. |fNSigTPC|<3, as an independent purity cut.\n");

    // ---------------------------------------------------------------
    // STEP 8: purity check for the FINAL selection (ellipse + PID), not
    // the ellipse-only sample checked in Step 5. This is a first concrete
    // number for "how much of our clean sample is confirmed-material
    // conversion" -- necessary but not sufficient on its own, since real
    // conversions also happen at the beampipe (~3cm) and elsewhere in the
    // tracking volume, not only in the 15-35cm ITS window.
    // ---------------------------------------------------------------
    Long64_t nFinalAndRadius = tree->GetEntries(gammaPlusPid && gammaRadius);
    printf("\n--- Purity cross-check on FINAL selection (ellipse + PID) ---\n");
    printf("Final selected sample: %lld\n", nGammaPlusPid);
    printf("Of those, %lld (%.1f%%) also fall in the 15-35 cm radius window.\n",
           nFinalAndRadius, 100.0*nFinalAndRadius/nGammaPlusPid);
    printf("This is a lower bound on confirmed-material-conversion purity,\n");
    printf("not a full contamination fraction -- real conversions also occur\n");
    printf("at the beampipe and elsewhere, not only in this one radius window.\n");
}