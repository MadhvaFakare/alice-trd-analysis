// Sum entries across every DF_*/trdpid in the writer output.
// Usage: root -l 'count_trdpid.C'

void count_trdpid(const char *inFile =
    "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root") {

    TFile *f = TFile::Open(inFile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", inFile); return; }

    Long64_t total = 0;
    int nDF = 0, nNonEmpty = 0;

    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;
        ++nDF;

        TDirectory *dir = (TDirectory*)f->Get(name);
        if (!dir) continue;

        TTree *t = (TTree*)dir->Get("trdpid");
        if (!t) {
            printf("%s : no tree named trdpid\n", name.Data());
            continue;
        }
        Long64_t n = t->GetEntries();
        total += n;
        if (n > 0) ++nNonEmpty;
    }

    printf("\nDF folders found      : %d\n", nDF);
    printf("DF folders with rows  : %d\n", nNonEmpty);
    printf("Total entries         : %lld\n\n", total);

    if (total > 0) {
        printf("The tree has per-track data. Species breakdown:\n");
        // chain everything for a quick species count
        TChain ch("trdpid");
        TIter next2(f->GetListOfKeys());
        while ((key = (TKey*)next2())) {
            TString name = key->GetName();
            if (!name.BeginsWith("DF_")) continue;
            ch.Add(TString(inFile) + "/" + name + "/trdpid");
        }
        printf("  chained entries : %lld\n", ch.GetEntries());
        printf("  electrons (0)   : %lld\n", ch.GetEntries("fSpecies==0"));
        printf("  pions     (1)   : %lld\n", ch.GetEntries("fSpecies==1"));
        printf("  unknown  (-1)   : %lld\n", ch.GetEntries("fSpecies==-1"));
    } else {
        printf("No rows were written. The schema exists but the table is empty.\n");
    }
}