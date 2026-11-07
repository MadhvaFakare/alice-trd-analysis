// =============================================================================
//  trd_io.h  --  plumbing only. Opens trdpid_tree.root and reads one track at
//                a time. No physics in here, so you never need to edit it.
//
//  Use it at the top of your own macro (same folder):
//      #include "trd_io.h"
// =============================================================================
#ifndef TRD_IO_H
#define TRD_IO_H

#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TString.h>
#include <iostream>

const char* const kDefaultFile =
  "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root";

const int kNLayers  = 6;
// fSpecies codes from PIDStudy (-1 means untagged). Not called kElectron:
// ROOT already uses that name for the PDG code 11 (TPDGCode.h).
const int kSpeciesElectron = 0;
const int kSpeciesPion     = 1;

// One track. After c->GetEntry(i), these hold the values of track i.
struct Track {
  Float_t q[3][kNLayers];  // q[window][layer]; window 0,1,2 = Q0,Q1,Q2; -1 = missing
  Float_t pt;
  Int_t   species;
};

// Opens the file, chains every DF_* folder, connects the branches to t.
TChain* OpenTrdTree(Track& t, const char* fileName = kDefaultFile)
{
  TFile* f = TFile::Open(fileName);
  if (!f || f->IsZombie()) { std::cout << "cannot open " << fileName << "\n"; return nullptr; }

  TChain* c = new TChain("trdpid");
  TIter next(f->GetListOfKeys());
  while (TKey* k = (TKey*)next()) {
    TString n = k->GetName();
    if (n.BeginsWith("DF_")) c->Add(Form("%s/%s/trdpid", fileName, n.Data()));
  }
  f->Close();

  for (int w = 0; w < 3; ++w)
    for (int l = 0; l < kNLayers; ++l)
      c->SetBranchAddress(Form("fQ%d%d", w, l), &t.q[w][l]);
  c->SetBranchAddress("fPt",      &t.pt);
  c->SetBranchAddress("fSpecies", &t.species);

  std::cout << "opened " << c->GetEntries() << " tracks\n";
  return c;
}

#endif