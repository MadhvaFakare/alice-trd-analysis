void learn_gammatag1(){
    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    TTree *tree = (TTree*)f->Get("DF_2336518081137728/O2tpcskimv0wde");
    if (!tree){
        printf("Error, could not find the tree");
        return;
    }
    Long64_t nTotal = tree->GetEntries();
    printf("Total entries: %lld\n\n", nTotal);
    // Step 1: What PID hypothesis actually exists in the tree
    TCanvas *cPID = new TCanvas("cPID","fPidIndex",800,600);
    tree->Draw("fPidIndex>>hPID(20,-0.5,19.5)");
    TH1F *hPID = (TH1F*)gDirectory->Get("hPID");
    hPID->SetTitle("fPidIndex;PID Index; Counts");
    cPID->SaveAs("pidindex_check1.png");
    printf("--- fPidIndex populated bins ---\n");
    for (int i = 1; i <= hPID->GetNbinsX();i++) {
        double content = hPID->GetBinContent(i);
        if (content > 0) {
            printf("index %0.f: %0.f Entries\n", hPID->GetBinCenter(i), content);
        }
    }
    // Step 2: Kinematic gamma-candidate region, no mass involved. Made the selection (using TCut) of the 
    // parameters qTV0 and AlphaV0, by selecting the branches and then applied the formula
    // |fAlphaV0| < 0.95 && fQtV0 < 0.05*sqrt(1 - (fAlphaV0)² / 0.95² )) 
    // Refer to https://arxiv.org/pdf/1803.09857 pg-7
    double qTmax = 0.05; 
    double alphaMax = 0.95;
    TCut gammaKinematic = Form("TMath::Abs(fAlphaV0)<%f && fQtV0 < %f*TMath::Sqrt(1 - (fAlphaV0*fAlphaV0)/(%f*%f))", 
        alphaMax,qTmax,alphaMax,alphaMax);
    Long64_t nGammaKinematic = tree->GetEntries(gammaKinematic);
    printf("---Kinematic Gamma Candidate region (ALICE PCM ellipse: qTmax=%.2f, alphaMax=%.2f)---\n", qTmax, alphaMax);
    printf("Entries: %lld (Fraction of total: %.6f)\n\n",nGammaKinematic,(double)nGammaKinematic/nTotal);
    // step 3: Does fGammaPsiPair actually discriminates in the region?
    TCanvas *cPsi = new TCanvas("cPsi","fGammaPsiPair",800,600);
    tree->Draw("fGammaPsiPair>>hPsiAll(200,-2,2)");
    TH1F *hPsiAll = (TH1F*)gDirectory->Get("hPsiAll");
    hPsiAll->SetLineColor(kBlack);
    hPsiAll->SetTitle("fGammaPsiPair; Psi pair; Counts");
    tree->SetLineColor(kMagenta);
    tree->Draw("fGammaPsiPair>>hPsiKinGamma(200,-2,2)",gammaKinematic, "same");
    cPsi->SaveAs("PsiPair_Check.png");
    printf("---fGammaPsiPair check---\n");
    printf("Black = full sample. Magenta = kinematic Gamma-candidate subset only.\n");
    // Step 4: Radius Cross check, independently of PID and mass entirely
    TCanvas *cRad = new TCanvas("cRad","fRadiusV0 comparision",800,600);
    tree->SetLineColor(kBlack);
    tree->Draw("fRadiusV0>>hRadAll(200,0,200)");
    TH1F *hRadAll = (TH1F*)gDirectory->Get("hRadAll");
    hRadAll->SetTitle("fRadiusV0;Radius (cm); Counts");
    tree->SetLineColor(kMagenta);
    tree->Draw("fRadiusV0>>hRadGamma(200,0,200)",gammaKinematic, "same");
    TH1F *hRadGamma = (TH1F*)gDirectory->Get("hRadGamma");
    hRadGamma->SetLineColor(kMagenta);
    cRad->SaveAs("Radius_Check.png");
    // Step 5: Take the corrected alpha and the qT back from the selection to the Armenteros-Podolanski plot
    TCanvas *cArmRefined = new TCanvas("cArmRefined","Armenteros with ALICE PCM gamma selection", 1000, 800);
    tree->Draw("fQtV0:fAlphaV0>>hAllRefined(200,-1,1, 200, 0,0.3)","","colz");
    TH2F *hAllRefined = (TH2F*)gDirectory->Get("hAllRefined");
    hAllRefined->SetTitle("Armenteros-Podolanski with ALICE PCM gamma selection;#alpha;q_{T} (GeV/c)");
    tree->SetMarkerColor(kMagenta);
    tree->SetMarkerStyle(20);
    tree->SetMarkerSize(0.4);
    tree->Draw("fQtV0:fAlphaV0", gammaKinematic, "same");
    cArmRefined->SaveAs("armenteros_refined_gammatag1.png");
    // optional-cross check against the material-radius clustering found
    TCut gammaRadius = "fRadiusV0>15 && fRadiusV0<35";
    Long64_t nBoth = tree->GetEntries(gammaKinematic && gammaRadius);
    printf("---Interpretation---\n");
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
    // Step 6: Test's the fix if fPIDIndex removes the contamination from the Armenteros-Podolanski
    TCut gammaPlusPid = gammaKinematic && "fPidIndex==0";
    Long64_t nGammaPlusPid = tree->GetEntries(gammaPlusPid);
    printf("\n --- Ellipse + fPidIndex==0---");
    printf("Entries: %lld (was %lld before the extra cut)\n\n",nGammaPlusPid,nGammaKinematic);
    TCanvas *cArmPid = new TCanvas("cArmPid","Armenteros with ellipse + PID cut", 1000,800);
    tree->Draw("fQtV0:fAlphaV0>>hAllPid(200,-1,1,200,0,0.3)","","colz");
    TH2F *hAllPid = (TH2F*)gDirectory->Get("hAllPid");
    hAllPid->SetTitle("Armenteros-Podolanski, ellipse + fPidIndex==0;#alpha;q_{T} (GeV/c)");
    tree->SetMarkerColor(kMagenta);
    tree->SetMarkerStyle(20);
    tree->SetMarkerSize(0.4);
    tree->Draw("fQtV0:fAlphaV0", gammaPlusPid, "same");
    cArmPid->SaveAs("armenteros-elipse-plus-pid.png");
    printf("Check: do the two blocks sitting on the horn\n");
    printf("If they dissaper, fPidIndex == 0 is doing real work\n");
    // Step 7: The fNSigTPC audit
    TCanvas *cNSig = new TCanvas("cNSig", "fNSigTPC by fPidIndex group", 800, 600);
    tree->SetLineColor(kRed);
    tree->Draw("fNSigTPC>>hNSig0(200,-10,10)", "fPidIndex==0");
    TH1F *hNSig0 = (TH1F*)gDirectory->Get("hNSig0");
    hNSig0->SetLineColor(kRed);
    hNSig0->SetTitle("fNSigTPC by fPidIndex group (red=0, blue=2, green=4);n#sigma_{TPC};Counts");
    tree->SetLineColor(kBlue);
    tree->Draw("fNSigTPC>>hNSig2(200,-10,10)","fPidIndex==2","same");
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
    // STEP 8: purity check for the FINAL selection (ellipse + PID), not
    // the ellipse-only sample checked in Step 5. This is a first concrete
    // number for "how much of our clean sample is confirmed-material
    // conversion" -- necessary but not sufficient on its own, since real
    // conversions also happen at the beampipe (~3cm) and elsewhere in the
    // tracking volume, not only in the 15-35cm ITS window.
    Long64_t nFinalAndRadius = tree->GetEntries(gammaPlusPid && gammaRadius);
    printf("\n--- Purity cross-check on FINAL selection (ellipse + PID) ---\n");
    printf("Final selected sample: %lld\n", nGammaPlusPid);
    printf("Of those, %lld (%.1f%%) also fall in the 15-35 cm radius window.\n",
           nFinalAndRadius, 100.0*nFinalAndRadius/nGammaPlusPid);
    printf("This is a lower bound on confirmed-material-conversion purity,\n");
    printf("not a full contamination fraction -- real conversions also occur\n");
    printf("at the beampipe and elsewhere, not only in this one radius window.\n");
}