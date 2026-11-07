const int kNLayer = 6;
const int kNSlice = 3;

void trd_likelihood01(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root"){
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) return;

    long nE = 0, nP = 0;

    TIter next(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())){
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;
        TDirectory* dir = (TDirectory*)key->ReadObj();
        TTree* t = (TTree*)dir->Get("trdpid");
        if (!t) continue;                                   // fix 5

        TLeaf* lSp = t->GetLeaf("fSpecies");
        TLeaf* lPt = t->GetLeaf("fPt");
        TLeaf* lQ[kNSlice][kNLayer];
        bool ok = (lSp && lPt);
        for (int s = 0; s < kNSlice && ok; s++){            // fix 6: && ok
            for (int l = 0; l < kNLayer; l++){
                lQ[s][l] = t->GetLeaf(Form("fQ%d%d", s, l));
                if (!lQ[s][l]) { ok = false; break; }
            }
        }                                                   // fix 6: closes the s loop
        if (!ok) continue;

        t->SetBranchStatus("*", 0);                         // read only what is needed
        t->SetBranchStatus("fSpecies", 1);
        t->SetBranchStatus("fPt", 1);
        for (int s = 0; s < kNSlice; s++){
            for (int l = 0; l < kNLayer; l++){
                t->SetBranchStatus(Form("fQ%d%d", s, l), 1);   // fix 4: 1, not l
            }
        }

        for (Long64_t i = 0; i < t->GetEntries(); i++){
            t->GetEntry(i);                                 // without this, GetValue() returns stale values
            int sp = (int)lSp->GetValue();
            if (sp != 0 && sp != 1) continue;               // fix 3
            if (sp == 0) nE++; else nP++;
        }
    }
    printf("Electrons %ld, pions %ld\n", nE, nP);
}