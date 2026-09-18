// List everything the PIDStudy task wrote, with class, dimensionality and
// entries. Run this before assuming Q0_El / Q0_Pi are all that is available —
// if a 2D Q0-vs-momentum histogram already exists, the momentum-binned study
// can be done without waiting for the flat tree.
//
// Usage:  root -l 'inspect_trd_output.C'
//         root -l 'inspect_trd_output.C("/path/to/AnalysisResults.root")'

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1.h"
#include "TTree.h"
#include "TClass.h"
#include "TString.h"

void walk(TDirectory *dir, const TString &path, int depth = 0) {
    TIter next(dir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TClass *cl = TClass::GetClass(key->GetClassName());
        const TString name = key->GetName();
        const TString full = path.Length() ? path + "/" + name : name;

        if (cl && cl->InheritsFrom(TDirectory::Class())) {
            printf("%*s[dir ] %s\n", depth * 2, "", name.Data());
            TDirectory *sub = (TDirectory*)key->ReadObj();
            walk(sub, full, depth + 1);
        }
        else if (cl && cl->InheritsFrom(TTree::Class())) {
            TTree *t = (TTree*)key->ReadObj();
            printf("%*s[TREE] %-32s entries=%lld branches=%d\n",
                   depth * 2, "", name.Data(), t->GetEntries(),
                   t->GetListOfBranches() ? t->GetListOfBranches()->GetEntries() : 0);
            if (t->GetListOfBranches()) {
                TIter nb(t->GetListOfBranches());
                TObject *b;
                while ((b = nb())) printf("%*s        branch: %s\n", depth * 2, "", b->GetName());
            }
        }
        else if (cl && cl->InheritsFrom(TH1::Class())) {
            TH1 *h = (TH1*)key->ReadObj();
            printf("%*s[%dD  ] %-32s entries=%.0f  bins=%d",
                   depth * 2, "", h->GetDimension(), name.Data(),
                   h->GetEntries(), h->GetNbinsX());
            if (h->GetDimension() > 1) printf(" x %d", h->GetNbinsY());
            printf("   x-axis: %s\n", h->GetXaxis()->GetTitle());
        }
        else {
            printf("%*s[%s] %s\n", depth * 2, "", key->GetClassName(), name.Data());
        }
    }
}

void inspect_trd_output(
    const char *inFile = "/home/madhvafakare/alice/PbPb2023pass5/trd/AnalysisResults.root") {
    TFile *f = TFile::Open(inFile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", inFile); return; }
    printf("\n=== Contents of %s ===\n", inFile);
    walk(f, "");
    printf("\nLook for: any 2D histogram with momentum on an axis, any TTree,\n");
    printf("and any Q1/Q2 or per-layer charge histograms beside Q0.\n");
    f->Close();
}