#include "TTree.h"
#include "Tcut.h"
void armenteros_by_pid1(){
    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/"
                            "LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    if (!f || f->IsZombie()){
        printf("Error opening file \n");
        return;
    } else {
        printf("File open's successfully!\n");
    }
    TTree *t = (TTree*)f->Get("DF_2336518081137728/O2tpcskimv0wde");
    if(!t){
        printf("Error opening tree\n");
        f->Close();
        return;
    } else {
        printf("Tree opens successfully!\n");
    }
    Long64_t nTotal = t->GetEntries();
    printf("Total entries: %lld\n\n", nTotal);
    // Selection uaing TCut
    TCut noCut = "";
    TCut pidElectron = "fPidIndex==0";
    TCut pidPion = "fPidIndex==2";
    TCut pidProton = "fPidIndex==4";
    
}