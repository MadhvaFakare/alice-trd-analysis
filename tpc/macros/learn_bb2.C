void learn_bb2(){
    TFile* f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root");
    if (!f || f->IsZombie()){
        cout<<"Error opeing file"<<endl;
        return;
    }
    TTree* tree = (TTree*)f->Get("DF_2336518081137728/O2tpcskimv0wde");
    if (!tree){
        cout<<"Error! No Tree found"<<endl;
        return;
    }
    /* o2::track::PID (ReconstructionDataFormats/PID.h)
    * Electron = 0, Muon = 1, Pion = 2, Kaon = 3, Proton = 4
    */
    const UChar_t kPidElectron = 0;
    const UChar_t kPidPion = 2;
    const UChar_t kPidProton = 4;
    TH2F* hAll = new TH2F("hAll","All species;#it{p} (GeV/#it{c}); TPC Signal(arb. units)", 150,0.05,5.0, 200, 0, 300);
    TH2F* hEl = new TH2F("hEl", "Electrons (V0: #gamma#rightarrowe^{+}e^{-});#it{p} (GeV/#it{c}); TPC Signal(arb. units)", 150,0.05,5.0, 200, 0, 300);
    TH2F* hPi = new TH2F("hPi", "Pions (V0: K^{0}_{s}#rightarrow#pi^{+}#pi^{-});#it{p} (GeV/#it{c}); TPC Signal(arb. units)", 150,0.05,5.0, 200, 0, 300);
    TH2F* hPr = new TH2F("hPr", "Proton (V0: #Lambda#rightarrowp#pi^{-});#it{p} (GeV/#it{c}); TPC Signal(arb. units)", 150,0.05,5.0, 200, 0, 300);    
    Float_t sig, mom;
    UChar_t pid;
    tree->SetBranchAddress("fTPCSignal", &sig);
    tree->SetBranchAddress("fTPCInnerParam", &mom);
    tree->SetBranchAddress("fPidIndex", &pid);
    for(Long64_t i = 0; i < tree->GetEntries();i++){
        tree->GetEntry(i);
        if (mom < 0.05 || mom > 5.0) continue;
        hAll->Fill(mom, sig);
        if(pid == kPidElectron) hEl->Fill(mom, sig);
        if(pid == kPidPion) hPi->Fill(mom, sig);
        if(pid == kPidProton) hPr->Fill(mom, sig);
    }
    TCanvas* c = new TCanvas("c","Species Seperation", 1600, 500);
    c->Divide(4,1); // 4 panels now
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    c->cd(1); gPad->SetLogx(); gPad->SetLogz(); gPad->SetRightMargin(0.15); hAll->Draw("COLZ");
    c->cd(2); gPad->SetLogx(); gPad->SetLogz(); gPad->SetRightMargin(0.15); hEl->Draw("COLZ");
    c->cd(3); gPad->SetLogx(); gPad->SetLogz(); gPad->SetRightMargin(0.15); hPi->Draw("COLZ");
    c->cd(4); gPad->SetLogx(); gPad->SetLogz(); gPad->SetRightMargin(0.15); hPr->Draw("COLZ");

    c->SaveAs("plots/step5_species.png");
    c->SaveAs("plots/step5_species.pdf");
    cout << "Done!" << endl;
}