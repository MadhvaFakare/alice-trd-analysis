static void collect(TTree* t, const char* cut, std::vector<double>& out){
    t->SetEstimate(t->GetEntries() + 1);
    Long64_t n = t->Draw("fQ00+fQ10+fQ20", cut, "goff");
    double* v = t->GetV1();
    for(Long64_t i = 0; i < n; i++) out.push_back(v[i]);
}
static void summary(const char* name, const std::vector<double>& q, double& mean, double& sd){
    mean = TMath::Mean(q.size(), q.data());
    sd = TMath::StdDev(q.size(), q.data());
    printf("   mean = %7.2f   std dev = %6.2f   error on mean = %5.2f\n",
         name, q.size(), mean, sd, sd / std::sqrt((double)q.size()));
}
void sumQ_layer0_1(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root", double pTMin = 0., double pTMax = 100.){
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()){
        printf("Error opening %s file\n",fname);
    } else {
        printf("File %s opens successfully!!\n",fname);
    }
    // layer 0 must have real measurements in all three slices
    // -1 is not an usable tracklet marker, each slides needs to be > 0
    TString base = Form("fQ00>0 && fQ10>0 && fQ20>0 && fPt>%g && fPt<%g", pTMin, pTMax);
    // Kangkan's species label: Electrons = 0, Pions = 1 and Unknown Particles = -1
    TString cutE = base + " && fSpecies==0";
    TString cutPi = base + " && fSpecies==1";
    // The file holds one trdpid tree per DF_* directory
    std::vector<double> qE, qPi;
    int nTrees = 0;
    TIter next(f->GetListOfKeys());
    while (TKey* key = (TKey*)next()){
        if (!TString(key->GetName()).BeginsWith("DF_")) continue;
        TTree* t = (TTree*)f->Get(Form("%s/trdpid", key->GetName()));
        if (!t) continue;
        collect(t, cutE.Data(), qE);
        collect(t, cutPi.Data(), qPi);
        ++nTrees;
    }
    if (qE.empty() || qPi.empty()) { printf("No electrons or no pions selected, check the file/cuts\n"); return; }

}