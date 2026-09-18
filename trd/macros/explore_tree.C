void explore_tree(const char *inFile ="/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root"){
    TFile *f = TFile::Open(inFile);
    if (!f || f->IsZombie()){
        printf("Error opening file %s\n", inFile);
        return;
    }
    // Step 1. The whole table in one tree
    TChain *t = new TChain("trdpid");
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())){
        TString name = key->GetName();
        if (name.BeginsWith("DF_")){
            t->Add(TString(inFile)+"/"+name+"/trdpid");
        } 
    }
    printf("\n=== Step 1. the whole table in one tree === \n");
    printf("Total tracks = %lld\n", t->GetEntries());
    // Step 2. Get the speices of the particles 
    const Long64_t nEl = t->GetEntries("fSpecies==0");
    const Long64_t nPi = t->GetEntries("fSpecies==1");
    const Long64_t nUnk = t->GetEntries("fSpecies==-1");
    printf("\n === Step 2. Species Breakdown ===\n");
    printf("Electrons: %lld\n",nEl);
    printf("Pions: %lld\n",nPi);
    printf("Unknown: %lld\n",nUnk);
    printf("So only %.1f%% of tracks carries a species label\n", 100.0 * (nEl+nPi)/t->GetEntries());
    // Step 3. The -1 Charge value
    printf("\n=== Step 3. why we need fQ00 > 0 ===\n");
    printf("electrons all row :%lld\n", nEl);
    printf("electrons with the usable layer 0 :%lld\n",t->GetEntries("fSpecies==0 && fQ00>0"));
    // Shows the damage completely, same damage same species, with and without the marker
    TH1F hBad("hBad", "",400,-10,200);
    TH1F hGood("hGood", "", 400, -10,200);
    t->Project("hBad","fQ00","fSpecies==0");
    t->Project("hGood", "fQ00", "fQ00>0 && fSpecies==0");
    printf("mean fQ00, -1s is included : %.2f <--wrong\n", hBad.GetMean());
    printf("mean fQ00, -1s is excluded : %.2f <--correct\n", hGood.GetMean());
    // Step 4. Look at some actual rows
    // Scan prints the table. Filtering on fSpecies>=0 skips the Unknown. That is 95% of the file
    printf("\n=== Step 4. ten tagged tracks ===\n");
    t->Scan("fPt:fSpecies:fQ00:fQ10:fQ20", "fSpecies>=0 && fQ00>0","");
    // Step 5. Project fills a histogram from the tree with a cut applied, then the histogram's mean is the mean of the selected rows.
    printf("\n=== Step 5. layer-0 means, usable layers only ===\n");
    printf("%-10s %12s %12s %12s %12s \n", "slice", "mean e", "N e", "mean pi", "N pi");
    const char *var[3] = {"fQ00", "fQ10", "fQ20"};
    const char *label[3] = {"Q0","Q1","Q2"};
    for (int i = 0; i < 3; i++){
        TH1F hE(Form("hE%d",i), "",200,0,200);
        TH1F hP(Form("hP%d",i), "",200,0,200);
        t->Project(hE.GetName(), var[i], Form("fSpecies==0 && %s>0", var[i]));
        t->Project(hP.GetName(), var[i], Form("fSpecies==1 && %s>0", var[i]));
        printf("%-10s %12.2f %12.0f %12.2f %12.0f\n",
            label[i], hE.GetMean(), hE.GetEntries(),
                      hP.GetMean(), hP.GetEntries());
    }
    // Step 6. Tagged tracks by momentum 
    printf("\n=== Step 6. tagged tracks by momentum ===\n");
    printf("%-16s %10s %10s\n", "pT window", "electrons", "pions");
    const double edge[6] = {0.0, 0.5, 1.0, 2.0, 4.0, 20.0};
    for (int i = 0; i < 5; i++){
        TString cE = Form("fSpecies==0 && fQ00 && fPt>%f && fPt<%f", edge[i], edge[i+1]);
        TString cP = Form("fSpecies==1 && fQ00 && fPt>%f && fPt<%f", edge[i], edge[i+1]);
        printf("%4.1f - %4.1f GeV/c %10lld %10lld\n",
            edge[i], edge[i+1], t->GetEntries(cE), t->GetEntries(cP));
    }
    // Step 7. The per-track sum of all three slices
    printf("\n=== STEP 7: per-track sum Q0+Q1+Q2, layer 0 ===\n");
    TH1F hSumE("hSumE", "", 300, 0, 300);
    TH1F hSumP("hSumP", "", 300, 0, 300);
    t->Project("hSumE", "fQ00+fQ10+fQ20", "fSpecies==0 && fQ00>0 && fQ10>0 && fQ20>0");
    t->Project("hSumP", "fQ00+fQ10+fQ20", "fSpecies==1 && fQ00>0 && fQ10>0 && fQ20>0");
    printf("electrons: mean %.2f  (N %.0f)\n", hSumE.GetMean(), hSumE.GetEntries());
    printf("pions    : mean %.2f  (N %.0f)\n", hSumP.GetMean(), hSumP.GetEntries());
 
    printf("\nDone. Every number above came from one line of code on the tree.\n");
    printf("The pattern is always the same: a variable, and a cut.\n\n");
}