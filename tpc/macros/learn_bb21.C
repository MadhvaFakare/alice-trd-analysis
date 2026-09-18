/* 
Author: Madhva Fakare
Ploting the histogram from the Data File DF_2336518081137728 
PID tags are find inside "/home/madhvafakare/alice/O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h"

line 88

class PID
{
 public:
  // particle identifiers, continuos starting from 0
  typedef pid_constants::ID ID;

  static constexpr ID Electron = 0;
  static constexpr ID Muon = 1;
  static constexpr ID Pion = 2;
  static constexpr ID Kaon = 3;
  static constexpr ID Proton = 4;
  static constexpr ID Deuteron = 5;
  static constexpr ID Triton = 6;
  static constexpr ID Helium3 = 7;
  static constexpr ID Alpha = 8;
  ....}
*/


/*
Aleph BetheBloch function to fit
*/
Double_t BetheBloch(Double_t x, Double_t par){
  Double_t p = x[0];
  Double_t mass = par[5];
  Double_t bg = p / mass;
  Double_t beta = bg / TMath::Sqrt(1.0 + bg * bg);
  if (beta <= 0 || beta >=1 ) return 0;
  Double_t betaP4 = TMath::Power(beta, par[3]);
  Double_t invbgP5 = TMath::Power(1.0 / bg, par[4]);
  Double_t logArg = (par[2] + invbgP5);
  if (logArg <= 0) return 0;
  Double_t val = (par[0]/betaP4)*(par[1] - betaP4 - TMath::Log(logArg));
  if (val <=0 ) return 0;
  return val;
}

// The main root macro 

void learn_bb21{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(kBird);
  gRandom->SetSeed(42);

/*

Information of the PID data is from the 88th line of 

/home/madhvafakare/alice/O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h

*/
  const UChar_t kPidElectron = 0; 
  const UChar_t kPidPion = 2;
  const UChar_t kPidProton = 4;

// Cache if exists 

const char* kCache = "plots/bb21_cache.root";

TH2F* hAll = nullptr;
TH2F* hEl = nullptr;
TH2F* hPi = nullptr;
TH2F* hPr = nullptr;

TFile* fCache = TFile::Open(kCache);
if (fCache && !fCache->IsZombie()){
  hAll = (TH2F*)fCache->Get('hAll'); if (hAll)->SetDirectory(0);
  hEl  = (TH2F*)fCache->Get('hEl');  if (hEl)->SetDirectory(0);
  hPi = (TH2F*)fCache->Get('hPi');  if (hPi)->SetDirectory(0);
  hPr  = (TH2F*)fCache->Get('hPr');  if (hPr)->SetDirectory(0);
  fCache->Close();
  if(hAll && hEl && hPi && hPr) {
    cout<<"Loaded histogram from cache."<<endl;
  } else {
     hAll = hEl = hPi = hPr = nullptr;
  }
}

if (!hAll) {
    hAll = new TH2F("hAll", "All species;#it{p} (GeV/#it{c});TPC Signal (arb. units)",150, 0.05, 5.0, 200, 0, 300);
    hEl  = new TH2F("hEl",  "Electrons (V0: #gamma#rightarrowe^{+}e^{-});#it{p} (GeV/#it{c});TPC Signal (arb. units)",150, 0.05, 5.0, 200, 0, 300);
    hPi  = new TH2F("hPi",  "Pions (V0: K^{0}_{s}#rightarrow#pi^{+}#pi^{-});#it{p} (GeV/#it{c});TPC Signal (arb. units)",150, 0.05, 5.0, 200, 0, 300);
    hPr  = new TH2F("hPr",  "Proton (V0: #Lambda#rightarrowp#pi^{-});#it{p} (GeV/#it{c});TPC Signal (arb. units)",150, 0.05, 5.0, 200, 0, 300);

    vector<TString> files = {
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzf_apass5_250514/AO2D_merge_LHC23zzf.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzg_apass5_250514/AO2D_merge_LHC23zzg.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzh_apass5_250514/AO2D_merge_LHC23zzh.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzi_apass5_250514/AO2D_merge_LHC23zzi.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzk_apass5_250514/AO2D_merge_LHC23zzk.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzl_apass5_250514/AO2D_merge_LHC23zzl.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzm_apass5_250514/AO2D_merge_LHC23zzm.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzn_apass5_250514/AO2D_merge_LHC23zzn.root",
            "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5/LHC23zzo_apass5_250514/AO2D_merge_LHC23zzo.root"
    };

    for(int iFile = 0; iFile < (int)files.size(); iFile++){
      cout << "File " << iFile+1 << "/9..." << endl;
            TFile* f = TFile::Open(files[iFile]);
            if (!f || f->IsZombie()) { cout << "Skipping." << endl; continue;
    }
    

}



}


