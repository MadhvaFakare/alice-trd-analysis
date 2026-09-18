// learn_armenteros_verify.C
//
// Purpose: verify which V0 species (Lambda, AntiLambda, K0s, gamma) corresponds
// to which peak in the Armenteros-Podolanski plot, using fMass as a cross-check
// against the alpha/qT kinematics.
//
// Step 1: inspect fMass on its own, to see what hypothesis it was computed under.
// Step 2: overlay the Armenteros plot with color-coded mass-window selections.
//
// Run with: root -l learn_armenteros_verify.C

void learn_armenteros_verify() {

    // ---- path confirmed from ROOT Object Browser ----
    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    TTree *tree = (TTree*)f->Get("DF_2336518081137728/O2tpcskimv0wde");

    if (!tree) {
        printf("ERROR: could not find tree, check the path in gDirectory\n");
        return;
    }

    // ---------------------------------------------------------------
    // STEP 1: inspect fMass alone first, before assuming any hypothesis
    // ---------------------------------------------------------------
    TCanvas *cMass = new TCanvas("cMass", "fMass distribution", 900, 600);
    tree->Draw("fMass>>hMass(200,0,1.3)");
    TH1F *hMass = (TH1F*)gDirectory->Get("hMass");
    hMass->SetTitle("fMass;Mass (GeV/c^{2});Counts");
    cMass->SaveAs("fMass_check.png");

    printf("\n--- fMass check ---\n");
    printf("If you see peaks near 0.0 (gamma), 0.498 (K0s), and 1.116 (Lambda),\n");
    printf("fMass is a generic reconstructed mass and species-specific cuts below are valid.\n");
    printf("If you see only ONE peak, fMass is a single-hypothesis mass\n");
    printf("(commonly the gamma/electron-pair hypothesis for this kind of skim)\n");
    printf("and the Lambda/K0s cuts below will NOT be meaningful as written.\n\n");

    // ---------------------------------------------------------------
    // STEP 2: Armenteros plot with color-coded candidate mass windows
    // ---------------------------------------------------------------
    TCanvas *cArm = new TCanvas("cArm", "Armenteros-Podolanski (mass-tagged)", 1000, 800);

    // full population in gray as background reference
    tree->Draw("fQtV0:fAlphaV0>>hAll(200,-1,1,200,0,0.3)", "", "colz");
    TH2F *hAll = (TH2F*)gDirectory->Get("hAll");
    hAll->SetTitle("Armenteros-Podolanski;#alpha;q_{T} (GeV/c)");

    // ---------------------------------------------------------------
    // DIAGNOSTIC: is the white space between the arc/horns truly empty,
    // or is a faint population being washed out by the linear color scale?
    // ---------------------------------------------------------------
    Long64_t nTotal = tree->GetEntries();
    Long64_t nWhitePatch = tree->GetEntries("abs(fAlphaV0)<0.2 && fQtV0>0.12 && fQtV0<0.18");
    printf("\n--- White-space diagnostic ---\n");
    printf("Entries in patch |alpha|<0.2, 0.12<qT<0.18 (looks white in linear colz): %lld\n", nWhitePatch);
    printf("Total entries: %lld  (fraction: %.6f)\n", nTotal, (double)nWhitePatch/nTotal);
    printf("If this is 0 or tiny, the white patch is genuinely empty.\n");
    printf("If it is a meaningful number, logz below should reveal it.\n\n");

    // Candidate mass windows -- ONLY trust these if Step 1 showed multiple peaks
    TCut cutLambda   = "TMath::Abs(fMass-1.115683)<0.006";  // Lambda -> p pi-
    TCut cutAntiLam  = cutLambda; // same mass, distinguished by alpha sign below
    TCut cutK0s      = "TMath::Abs(fMass-0.497611)<0.010";  // K0s -> pi+ pi-
    TCut cutGamma    = "TMath::Abs(fMass)<0.020";           // gamma conversion

    tree->SetMarkerColor(kRed);
    tree->Draw("fQtV0:fAlphaV0", cutLambda && "fAlphaV0>0", "same");   // Lambda: alpha>0
    tree->SetMarkerColor(kBlue);
    tree->Draw("fQtV0:fAlphaV0", cutLambda && "fAlphaV0<0", "same");   // AntiLambda: alpha<0
    tree->SetMarkerColor(kGreen+2);
    tree->Draw("fQtV0:fAlphaV0", cutK0s, "same");
    tree->SetMarkerColor(kMagenta);
    tree->Draw("fQtV0:fAlphaV0", cutGamma, "same");

    cArm->SaveAs("armenteros_masstagged_linear.png");

    // redraw the SAME canvas with a log color scale to check whether a faint
    // low-level population is hiding in what looked like white/empty space above
    cArm->SetLogz(1);
    cArm->Modified();
    cArm->Update();
    cArm->SaveAs("armenteros_masstagged_logz.png");

    printf("Red    = Lambda candidates (alpha > 0)\n");
    printf("Blue   = AntiLambda candidates (alpha < 0)\n");
    printf("Green  = K0s candidates\n");
    printf("Magenta= gamma candidates\n");
    printf("\nCheck: does red land on your right-hand peak (~+0.65)?\n");
    printf("Does blue land on your left-hand peak (~-0.65)?\n");
    printf("Does green trace an arc through the central plateau?\n");
    printf("Does magenta cluster tightly near alpha=0 at LOW qT?\n");
}