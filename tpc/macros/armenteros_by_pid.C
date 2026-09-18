void armenteros_by_pid() {
    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/"
                           "LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    if (!f || f->IsZombie()) {
        printf("Error opening file\n");
        return;
    }
    TTree *t = (TTree*)f->Get("DF_2336518081137728/O2tpcskimv0wde");
    if (!t) {
        printf("Error: could not find tree\n");
        return;
    }
    Long64_t nTotal = t->GetEntries();
    printf("Total entries: %lld\n\n", nTotal);
    // ---- selections as TCut, same style as learn_gammatag1.C ----
    TCut noCut       = "";
    TCut pidElectron = "fPidIndex==0";
    TCut pidPion     = "fPidIndex==2";
    TCut pidProton   = "fPidIndex==4";
    TCut cuts[4]   = { noCut, pidElectron, pidPion, pidProton };
    const char *names[4]  = { "hAll", "hP0", "hP2", "hP4" };
    const char *titles[4] = { "No selection",
                              "fPidIndex==0 (electron)",
                              "fPidIndex==2 (pion)",
                              "fPidIndex==4 (proton)" };
    TCanvas *c = new TCanvas("cArmPid", "Armenteros by PID", 1200, 900);
    c->Divide(2,2);
    for (int i = 0; i < 4; i++) {
        c->cd(i+1);
        gPad->SetLogz();
        gPad->SetRightMargin(0.15);
        t->Draw(Form("fQtV0:fAlphaV0>>%s(200,-1,1,200,0,0.3)", names[i]),
                cuts[i], "colz");
        TH2F *h = (TH2F*)gDirectory->Get(names[i]);
        if (!h) {
             printf("could not retrieve %s\n", names[i]);
              continue; 
            }
        h->SetTitle(Form("%s;#alpha;q_{T} (GeV/c)", titles[i]));
        TString cs = cuts[i].GetTitle();
        Long64_t n = cs.Length() ? t->GetEntries(cs) : nTotal;
        printf("%-28s : %12lld entries (%.2f%% of total)\n",
               titles[i], n, 100.0*n/nTotal);
        gPad->Modified();
        gPad->Update();
    }
    c->SaveAs("armenteros_by_pidindex.png");
}