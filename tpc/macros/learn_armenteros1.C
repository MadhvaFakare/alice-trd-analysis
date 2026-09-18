void learn_armenteros1(){
    //---------------- Style------------------- 
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gRandom->SetSeed(42);
    //-----------------PID Index---------------
    const uChar_t kPidElectron = 0;
    //-----------------Histograms--------------
    // hArm: all V0 tracks --full Armenteros plot
    // hArmEl: electrons only -- where fPidIndex = 0 sits
    TH2F* hArm = new TH2F("hArm","Armenteros-Podolanski - All V0 tracks;"
                        "#alpha_{V0} = (P_{L}^{+} - P_{L}^{-})/(P_{L}^{+} + P_{L}^{-})"
                        "#it{q}_{T} (GeV/#it{c}",
                        200, -1.0, 1,0, 200, 0.0, 0.25);
    TH2F* hArmEL = new TH2F("hArmEl","Armenteros-Podolanski - V0 electrons (fPidIndex=0);"
                        "#alpha_{V0} = (P_{L}^{+} - P_{L}^{-})/(P_{L}^{+} + P_{L}^{-})"
                        "#it{q}_{T} (GeV/#it{c}",
                        200, -1.0, 1,0, 200, 0.0, 0.25);
    //-----------------Open File----------------
    TFile* f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/"
                            "LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    if (!f || f->IsZombie()){
        cout<<"Error: could not open the file."<<endl;
        return;
    }
    cout<<"File opened successfully";
    //-----------------Loop over DF_folders-----
    TIter next(f->GetListOfKeys());
    TKey* key;
    int nFolders = 0;

    while ((key = (TKey*)next())){
        TString name = key->GetName();
        if (!name.BeginWith("DF_")) continue;

        TTree* tree = (TTree*)f->Get(name + "O2tpcskimv0wde");
        if(!tree) continue;
        nFolders++;
        //Only read the branch that we need
        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("fAlphaV0")
    }
}