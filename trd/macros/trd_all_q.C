void trd_all_q(const char* inFile = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root"){
    TFile *f = TFile::Open(inFile);
    if (!f || f->IsZombie()) {printf("Error opening File\n"); return;} else {printf("Successfuly open file\n");}

}