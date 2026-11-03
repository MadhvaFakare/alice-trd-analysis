void trd_separation01(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root"){
   TFile *f = TFile::Open(fname); // Open
   if (!f || f->IsZombie()) return;
   TH1D* hE = new TH1D("hE", "Electrons Q0 layer 3;charge;counts",64,0,128);
   TH1D* hP = new TH1D("hP", "Pions Q0 layer 3;charge;counts",64,0,128);
    TIter next(f->GetListOfKeys()); // List
    TKey* Key;
    while((Key = (TKey*)next())){ 
        TString name = Key->GetName(); // Filter 
        if (!name.BeginsWith("DF_")) continue; // skip and continue if there are no names that begins wih DF_
        TDirectory* dir = (TDirectory*)Key->ReadObj(); 
        TTree* t = (TTree*)dir->Get("trdpid");
        if (!t) continue;
        TLeaf* lSp = t->GetLeaf("fSpecies");
        TLeaf* lQ = t->GetLeaf("fQ03"); // Slice Q0 layer 3
        if (!lSp || !lQ) continue;
        Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; i++){
            t->GetEntry(i); // load track i into leaves
            int sp = (int)lSp->GetValue(); // 0 = electrons, 1 = pions
            double q = lQ->GetValue(); // charge of the track in layer 3
            if (q<=0) continue;
            if (sp==0) hE->Fill(q);
            if (sp==1) hP->Fill(q);
        }
    }
    printf("entries: electrons %.0f, pions %.0f\n", hE->GetEntries(), hP->GetEntries());
    // normalize: each histogram sums to 1, so shapes are comparable
    hE->Scale(1.0 / hE->Integral());
    hP->Scale(1.0 / hP->Integral());

    hE->SetLineColor(kRed + 1);
    hP->SetLineColor(kBlue + 1);

    // make the y axis tall enough for whichever histogram is higher
    double m = hE->GetMaximum();
    if (hP->GetMaximum() > m) m = hP->GetMaximum();
    hE->SetMaximum(1.2 * m);

    hE->Draw("hist");
    hP->Draw("hist same");
    
}