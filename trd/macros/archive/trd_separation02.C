void trd_separation02(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root"){
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) return;
    gStyle->SetOptStat("0");
    TH1D* hE = new TH1D("hE", "Electrons Q0 in layer 3;charge;counts",64,0,128);
    TH1D* hP = new TH1D("hP", "Pions Q0 in layer 3;charge;counts",64,0,128);
    TIter next(f->GetListOfKeys());
    TKey* key;
    while((key = (TKey*)next())){
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;
        TDirectory* dir = (TDirectory*)key->ReadObj();
        TTree* t = (TTree*)dir->Get("trdpid");
        if (!t) continue;
        TLeaf* lSp = t->GetLeaf("fSpecies"); // Species looks at all the particles including Electrons, Pions and Unknown particles 
        TLeaf* lQ = t->GetLeaf("fQ03"); // Specifically targeting Q03 for all the species particles 
        if (!lSp || !lQ) continue;
        Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; i++){
            t->GetEntry(i);
            int sp = (int)lSp->GetValue();
            double q = lQ->GetValue();
            if (q <= 0) continue; // remove the unkown particles Kangkan's pid file has -1 as UNknown particles
            if (sp == 0) hE->Fill(q); // Electrons are tagged as 0
            if (sp == 1) hP->Fill(q); // Pions are tagged as 1
        }
    }
    printf("Entries: electrons %.0f and pions %.0f\n", hE->GetEntries(), hP->GetEntries());
    // normalize both the histogram
    hE->Scale(1.0/hE->Integral());
    hP->Scale(1.0/hP->Integral());
    hE->SetLineColor(kRed + 1);
    hP->SetLineColor(kBlue + 1);
    // make the y-axis tall enough for whichever histogram is higher
    double m = hE->GetMaximum();
    if (hP->GetMaximum() > m) m = hP->GetMaximum();
    hE->SetMaximum(1.2 * m);
    hE->Draw("hist"); hP->Draw("hist same");
    TLegend* leg = new TLegend(0.65,0.72,0.88,0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hE,"Electron","l");
    leg->AddEntry(hP,"Pion","l");
    leg->Draw();
}