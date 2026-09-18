// purify_electrons.C
//
// Purpose: apply the validated photon-conversion electron selection
// (ALICE PCM ellipse in (alpha, qT) + fPidIndex==0) across ALL runs and
// ALL DF_* folders, and write the surviving candidates into a single
// output file, purified_electrons.root, keeping every original branch.
//
// Selection being applied (validated in learn_gammatag.C, see the
// session-summary PDF, Sections 8-11):
//   1) qT < 0.05 * sqrt(1 - alpha^2 / 0.95^2)   [ALICE PCM ellipse]
//   2) fPidIndex == 0                            [electron hypothesis]
//
// Known caveat carried over from the summary, Section 12: whether
// fPidIndex assignment is conditioned on an upstream preselection is
// still an open question for Guido/Anton/Kangkan. This script applies
// the selection as validated so far; revisit if their answer changes
// the picture.
//
// Usage:
//   root -l -q macros/purify_electrons.C
// or with a different base directory:
//   root -l -q 'macros/purify_electrons.C("/some/other/path")'

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include <vector>

// -------------------------------------------------------------------
// Find every DF_*/O2tpcskimv0wde inside one AO2D file and return the
// full chain-addressable paths ("file.root/DF_xxx/O2tpcskimv0wde").
// -------------------------------------------------------------------
std::vector<TString> findV0Trees(const TString &filePath) {
    std::vector<TString> paths;
    TFile *f = TFile::Open(filePath);
    if (!f || f->IsZombie()) {
        printf("  WARNING: could not open %s, skipping\n", filePath.Data());
        if (f) f->Close();
        return paths;
    }
    TIter nextKey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextKey())) {
        TString kname = key->GetName();
        if (!kname.BeginsWith("DF_")) continue;
        TDirectory *dir = (TDirectory*)f->Get(kname);
        if (!dir) continue;
        if (dir->GetListOfKeys()->FindObject("O2tpcskimv0wde")) {
            paths.push_back(filePath + "/" + kname + "/O2tpcskimv0wde");
        }
    }
    f->Close();
    return paths;
}

void purify_electrons(const char *baseDir =
        "/home/madhvafakare/alice/PbPb2023pass5/PbPb2023pass5") {

    // ---------------------------------------------------------------
    // STEP 1: discover AO2D_merge_*.root files in run subdirectories
    // ---------------------------------------------------------------
    std::vector<TString> files;
    TSystemDirectory base("base", baseDir);
    TList *entries = base.GetListOfFiles();
    if (!entries) {
        printf("ERROR: cannot read directory %s\n", baseDir);
        return;
    }
    TIter nextEntry(entries);
    TSystemFile *entry;
    while ((entry = (TSystemFile*)nextEntry())) {
        TString name = entry->GetName();
        if (name == "." || name == "..") continue;
        if (!entry->IsDirectory()) continue;
        // Skip test/subset merges: LHC23zzo_apass5_250514_fewRuns contains
        // a subset of the runs already in LHC23zzo_apass5_250514, so
        // including both double-counts those entries. Verify with the
        // fRunNumber overlap check if unsure.
        if (name.EndsWith("_fewRuns")) {
            printf("Skipping subset directory: %s\n", name.Data());
            continue;
        }
        // look inside each run folder for AO2D_merge_*.root
        TString runDir = TString(baseDir) + "/" + name;
        TSystemDirectory sub("sub", runDir);
        TList *subEntries = sub.GetListOfFiles();
        if (!subEntries) continue;
        TIter nextSub(subEntries);
        TSystemFile *subEntry;
        while ((subEntry = (TSystemFile*)nextSub())) {
            TString fname = subEntry->GetName();
            if (fname.BeginsWith("AO2D_merge_") && fname.EndsWith(".root")) {
                files.push_back(runDir + "/" + fname);
            }
        }
    }

    if (files.empty()) {
        printf("ERROR: no AO2D_merge_*.root files found under %s\n", baseDir);
        printf("Check the base directory path.\n");
        return;
    }

    printf("Found %zu AO2D file(s):\n", files.size());
    for (auto &fp : files) printf("  %s\n", fp.Data());
    printf("\n");

    // ---------------------------------------------------------------
    // STEP 2: build one chain over every DF_*/O2tpcskimv0wde found
    // ---------------------------------------------------------------
    TChain *chain = new TChain("O2tpcskimv0wde");
    int nTreesTotal = 0;
    for (auto &fp : files) {
        std::vector<TString> treePaths = findV0Trees(fp);
        printf("%s : %zu DF folder(s) with O2tpcskimv0wde\n",
               gSystem->BaseName(fp.Data()), treePaths.size());
        for (auto &tp : treePaths) {
            chain->Add(tp);
            nTreesTotal++;
        }
    }
    printf("\nChained %d tree(s) total.\n", nTreesTotal);

    Long64_t nTotal = chain->GetEntries();
    printf("Total entries across all runs: %lld\n\n", nTotal);
    if (nTotal == 0) {
        printf("ERROR: chain is empty, aborting.\n");
        return;
    }

    // ---------------------------------------------------------------
    // STEP 3: the validated selection
    // ---------------------------------------------------------------
    double qTmax = 0.05;
    double alphaMax = 0.95;
    TCut ellipse = Form(
        "TMath::Abs(fAlphaV0)<%f && fQtV0 < %f*TMath::Sqrt(1 - (fAlphaV0*fAlphaV0)/(%f*%f))",
        alphaMax, qTmax, alphaMax, alphaMax);
    TCut electronPid = "fPidIndex==0";
    TCut purified = ellipse && electronPid;

    // ---------------------------------------------------------------
    // STEP 4: per-run bookkeeping BEFORE the merge, so the thesis has
    // run-by-run numbers, not only one global figure
    // ---------------------------------------------------------------
    printf("--- Per-file bookkeeping ---\n");
    printf("%-40s %12s %12s %8s\n", "file", "entries", "selected", "frac");
    for (auto &fp : files) {
        std::vector<TString> treePaths = findV0Trees(fp);
        Long64_t nFile = 0, nSel = 0;
        for (auto &tp : treePaths) {
            TChain tmp("O2tpcskimv0wde");
            tmp.Add(tp);
            nFile += tmp.GetEntries();
            nSel  += tmp.GetEntries(purified);
        }
        printf("%-40s %12lld %12lld %7.2f%%\n",
               gSystem->BaseName(fp.Data()), nFile, nSel,
               nFile > 0 ? 100.0 * nSel / nFile : 0.0);
    }
    printf("\n");

    // ---------------------------------------------------------------
    // STEP 5: write the purified tree, all branches preserved
    // ---------------------------------------------------------------
    TFile *out = TFile::Open("purified_electrons.root", "RECREATE");
    if (!out || out->IsZombie()) {
        printf("ERROR: cannot create purified_electrons.root\n");
        return;
    }
    printf("Copying selected entries (this can take a while on the full set)...\n");
    TTree *pureTree = chain->CopyTree(purified);
    pureTree->SetName("purifiedElectrons");
    pureTree->SetTitle("Photon-conversion electron candidates, PCM ellipse + fPidIndex==0");

    Long64_t nPure = pureTree->GetEntries();
    pureTree->Write();

    // ---------------------------------------------------------------
    // STEP 6: summary + the radius purity indicator for the merged set
    // ---------------------------------------------------------------
    Long64_t nRadius = pureTree->GetEntries("fRadiusV0>15 && fRadiusV0<35");
    printf("\n--- Summary ---\n");
    printf("Total entries (all runs):        %lld\n", nTotal);
    printf("Purified electron candidates:    %lld  (%.2f%% of total)\n",
           nPure, 100.0 * nPure / nTotal);
    printf("In 15-35 cm material window:     %lld  (%.1f%% of purified)\n",
           nRadius, nPure > 0 ? 100.0 * nRadius / nPure : 0.0);
    printf("Reminder: the radius fraction is a LOWER bound on material-\n");
    printf("conversion purity, the window excludes the beampipe and inner\n");
    printf("ITS layers where real conversions also occur.\n");
    printf("\nOutput written: purified_electrons.root  (tree: purifiedElectrons)\n");

    out->Close();
}