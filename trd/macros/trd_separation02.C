void trd_separation02(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root"){
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) return;
    else printf("File name '%s' opens successfully \n",fname);
    TH1D* hE = new TH1D("hE","Electrons Q0 layer 3;charge;counts",64,0,128);
    TH1D* hP = new TH1D("hP", "Pions Q0 layer 3;charge;counts",64,0,128); 
    TIter next(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())){
        TString name = key->GetName();
        if(!name.BeginsWith("DF_")) continue;
        TDirectory* dir = (TDirectory*)key->ReadObj();
        TTree *t = (TTree*)dir->Get("trdpid");
        if (!t) continue;
        TLeaf* lsp = t->GetLeaf("fSpecies");
        TLeaf* lQ = t->GetLeaf("fQ03");
        if (!lsp || !lQ) continue;
        Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; i++){
            t->GetEntry(i);
            int sp = (int)lsp->GetValue();
            double q = lQ->GetValue();
            if (q <= 0) continue;
            if (sp == 0) hE->Fill(q);
            if (sp == 1) hP->Fill(q);
        }
    }
    printf("entries: electrons %0.f, pions %0.f \n", hE->GetEntries(), hP->GetEntries());
    hE->Scale(1.0/hE->Integral());
    hP->Scale(1.0/hP->Integral());
    hE->SetLineColor(kRed + 1);
    hP->SetLineColor(kBlue + 1);
    // make the y-axis taller for whichever histogram is higher
    double m = hE->GetMaximum();
    if (hP->GetMaximum() > m) m = hP->GetMaximum();
    hE->SetMaximum(1.2*m);
    hE->Draw("hist");
    hP->Draw("hist same");
    gPad->SaveAs("plots/trd_separation02_q0_layer3.png");
}