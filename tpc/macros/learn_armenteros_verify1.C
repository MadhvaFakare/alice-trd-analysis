void learn_armenteros_verify1(){
    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    if (!f || f->IsZombie()){
        cout<<"Error opening file"<<endl;
    } else {
        cout<<"File Open successfully"<<endl;
    }
    TTree *tree = (TTree*)f->Get("DF_2336518081137728/O2tpcskimv0wde");
    if (!tree) {
        cout<<"Error opening tree"<<endl;
    } else {
        cout<<"Tree opened successfully"<<endl;
    }
    // step 1. Inspect fMass alone first
    TCanvas *cMass = new TCanvas("cMass","fMass distribution", 900, 600);
    tree->Draw("fMass>>hMass(200,0,1.5)");
    TH1F *hMass = (TH1F*)gDirectory->Get("hMass");
    hMass->SetTitle("fMass;Mass (GeV/c^{2});Counts");
    cMass->SaveAs("fMass_check1.png");
    printf("Entries with fMass<0.02: %lld\n", tree->GetEntries("fMass<0.02"));
    // Step 2. Armenteros Plot with color-coded candidate mass window
    TCanvas *cArm = new TCanvas("cArm", "Armenteros-Podolanski", 1000, 800);
    // full population in grey as background reference 
    tree->Draw("fQtV0:fAlphaV0>>hAll(200,-1,1,200,0,0.3)","","colz");
    TH2F *hAll = (TH2F*)gDirectory->Get("hAll");
    hAll->SetTitle("Armenteros-Podolanski Plot#alpha;q_{T} (GeV/C)");
    //Verify if those white space are truely emplty
    Long64_t nTotal = tree->GetEntries();
    Long64_t nWhitePatch = tree->GetEntries("abs(fAlphaV0)<0.2 && fQtV0>0.12 && fQtV0<0.18");
    printf("\n--- White-space diagnostic ---\n");
    printf("Entries in patch |alpha|<0.2, 0.12<qT<0.18 (looks white in linear colz): %lld\n", nWhitePatch);
    printf("Total entries: %lld  (fraction: %.6f)\n", nTotal, (double)nWhitePatch/nTotal);
    printf("If this is 0 or tiny, the white patch is genuinely empty.\n");
    printf("If it is a meaningful number, logz below should reveal it.\n\n");
}