// Step 1: open the data and count electrons and pions in 1 < pT < 2 GeV/c.
// Run:  root -l -q lq_step1.C

#include "trd_io.h"

void lq_step1()
{
  Track t;
  TChain* c = OpenTrdTree(t);
  if (!c) return;

  int nEl = 0, nPi = 0;
  for (Long64_t i = 0; i < c->GetEntries(); ++i) {
    c->GetEntry(i);
    if (t.pt < 1.0 || t.pt >= 2.0) continue;
    if (t.species == kSpeciesElectron) ++nEl;
    if (t.species == kSpeciesPion)     ++nPi;
  }

  std::cout << "electrons: " << nEl << "   pions: " << nPi << "\n";
}