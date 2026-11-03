// trd_likelihood_felix.C
//
// LQ1D and LQ3D electron likelihood, following Felix Schlepper's thesis
// (Section 5.3, Eq. 5.3 to 5.7), applied to your flat tree.
//
// Per layer i, per slice j:   L_ij = P(Q_j|e) / ( P(Q_j|e) + P(Q_j|pi) )     (Eq. 5.4)
// LQ1D: one charge per layer, Q = Q0+Q1+Q2, factors over layers i             (Eq. 5.5)
// LQ3D: three charges per layer, factors over layers i and slices j           (Eq. 5.6, 5.7)
// Combined:  L = prod(L_e) / ( prod(L_e) + prod(1 - L_e) ) = 1 / (1 + exp(-S)),
//            with S = sum over factors of ln( P(Q|e) / P(Q|pi) ).
// The macro works with S (log space) so that nothing underflows.
//
// Also computed as a baseline: the mean charge per valid layer (like trd_separation04).
//
// Differences to the thesis (on purpose, because of ~150 electrons):
//   * pdfs are histograms with a pseudo-count (alpha), not fits. No low-charge cutoff.
//   * By default one pdf is shared by all layers (perLayer = false). Felix has one LUT per layer.
//   * One wide pT bin (all pT), not his four momentum intervals.
//   * 2-fold cross validation: tracks of DF_even are scored with pdfs built from DF_odd and
//     the other way round, so no track is scored by a pdf that contains itself.
//
// A layer counts as valid only if Q0, Q1 and Q2 are all > 0 (-1 is the "no tracklet" marker).
//
// Usage:
//   root -l 'macros/trd_likelihood_felix.C'
//   root -l 'macros/trd_likelihood_felix.C("/path/trdpid_tree.root", true)'   // per-layer pdfs

struct TrkL {
    int    sp;        // 0 electron, 1 pion
    int    fold;      // 0 or 1
    int    nv;        // number of valid layers
    bool   v[6];
    double q[6][3];
};

static int binOf(double q, double mx, int nb)
{
    int b = (int)(q / mx * nb);
    if (b < 0) b = 0;
    if (b >= nb) b = nb - 1;   // overflow goes into the last bin
    return b;
}

// pion efficiency at a given electron efficiency, using a score (higher = more electron-like)
static double pionEffAt(std::vector<double> e, const std::vector<double>& p, double effE)
{
    if (e.empty() || p.empty()) return -1;
    std::sort(e.begin(), e.end());
    size_t k = (size_t)std::floor((1.0 - effE) * e.size());
    if (k >= e.size()) k = e.size() - 1;
    double cut = e[k];
    long n = 0;
    for (double s : p) if (s > cut) n++;
    return (double)n / p.size();
}

void trd_likelihood_felix(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root",
                          bool perLayer = false, int minLayers = 4,
                          double ptMin = 0.0, double ptMax = 1e9)
{
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", fname); return; }
    gStyle->SetOptStat(0);

    const int nSlice = 3, nLayer = 6;
    const int    nB1 = 32;  const double max1 = 320.; // LQ1D: Q0+Q1+Q2, 10 units per bin
    const int    nB3 = 16;  const double max3 = 128.; // LQ3D: each slice, 8 units per bin
    const double alpha = 0.5;                          // pseudo-count per bin

    // ---------- 1. read all tagged tracks into memory ----------
    std::vector<TrkL> trk;
    int dfIndex = 0;
    TIter next(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;
        TDirectory* dir = (TDirectory*)key->ReadObj();
        TTree* t = (TTree*)dir->Get("trdpid");
        if (!t) continue;

        TLeaf* lSp = t->GetLeaf("fSpecies");
        TLeaf* lPt = t->GetLeaf("fPt");
        TLeaf* lQ[nSlice][nLayer];
        bool ok = (lSp && lPt);
        for (int s = 0; s < nSlice && ok; s++)
            for (int l = 0; l < nLayer; l++) {
                lQ[s][l] = t->GetLeaf(Form("fQ%d%d", s, l));
                if (!lQ[s][l]) { ok = false; break; }
            }
        if (!ok) continue;

        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("fSpecies", 1);
        t->SetBranchStatus("fPt", 1);
        for (int s = 0; s < nSlice; s++)
            for (int l = 0; l < nLayer; l++)
                t->SetBranchStatus(Form("fQ%d%d", s, l), 1);

        int fold = dfIndex % 2;
        dfIndex++;

        Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; i++) {
            t->GetEntry(i);
            int sp = (int)lSp->GetValue();
            if (sp != 0 && sp != 1) continue;
            double pt = lPt->GetValue();
            if (pt < ptMin || pt > ptMax) continue;

            TrkL tr;
            tr.sp = sp; tr.fold = fold; tr.nv = 0;
            for (int l = 0; l < nLayer; l++) {
                double q0 = lQ[0][l]->GetValue();
                double q1 = lQ[1][l]->GetValue();
                double q2 = lQ[2][l]->GetValue();
                tr.v[l] = (q0 > 0 && q1 > 0 && q2 > 0);
                tr.q[l][0] = q0; tr.q[l][1] = q1; tr.q[l][2] = q2;
                if (tr.v[l]) tr.nv++;
            }
            if (tr.nv < minLayers) continue;
            trk.push_back(tr);
        }
    }

    long cnt[2][2] = {{0, 0}, {0, 0}}; // [fold][species]
    for (const auto& tr : trk) cnt[tr.fold][tr.sp]++;
    printf("\nTracks used (>= %d valid layers): %zu\n", minLayers, trk.size());
    printf("  fold 0: %ld electrons, %ld pions\n", cnt[0][0], cnt[0][1]);
    printf("  fold 1: %ld electrons, %ld pions\n", cnt[1][0], cnt[1][1]);
    printf("  pdfs: %s, %d DF folders\n", perLayer ? "one per layer" : "shared by all layers", dfIndex);
    if (cnt[0][0] < 5 || cnt[1][0] < 5) { printf("Too few electrons in one fold, stopping.\n"); return; }

    // ---------- 2. build pdfs per fold: p[fold][species][layer][(slice)][bin] ----------
    static double p1[2][2][6][32];       // LQ1D
    static double p3[2][2][6][3][16];    // LQ3D
    for (int fo = 0; fo < 2; fo++) {
        double c1[2][6][32] = {};
        double c3[2][6][3][16] = {};
        double tot[2][6] = {};
        for (const auto& tr : trk) {
            if (tr.fold != fo) continue;
            for (int l = 0; l < nLayer; l++) {
                if (!tr.v[l]) continue;
                int li = perLayer ? l : 0;
                double sum = tr.q[l][0] + tr.q[l][1] + tr.q[l][2];
                c1[tr.sp][li][binOf(sum, max1, nB1)] += 1;
                for (int s = 0; s < nSlice; s++)
                    c3[tr.sp][li][s][binOf(tr.q[l][s], max3, nB3)] += 1;
                tot[tr.sp][li] += 1;
            }
        }
        for (int sp = 0; sp < 2; sp++)
            for (int li = 0; li < nLayer; li++) {
                for (int b = 0; b < nB1; b++)
                    p1[fo][sp][li][b] = (c1[sp][li][b] + alpha) / (tot[sp][li] + alpha * nB1);
                for (int s = 0; s < nSlice; s++)
                    for (int b = 0; b < nB3; b++)
                        p3[fo][sp][li][s][b] = (c3[sp][li][s][b] + alpha) / (tot[sp][li] + alpha * nB3);
            }
    }

    // ---------- 3. score every track with the pdfs of the other fold ----------
    // methods: 0 = LQ1D, 1 = LQ3D, 2 = mean charge (baseline)
    const char* mname[3] = {"LQ1D", "LQ3D", "mean charge"};
    std::vector<double> sE[3], sP[3];
    TH1D* hLe[3][2];
    for (int m = 0; m < 2; m++)
        for (int sp = 0; sp < 2; sp++)
            hLe[m][sp] = new TH1D(Form("hLe_%d_%d", m, sp),
                                  Form("%s;L_{e};fraction of tracks", mname[m]), 50, 0, 1);

    for (const auto& tr : trk) {
        int tf = 1 - tr.fold; // training fold
        double S1 = 0, S3 = 0, sumQ = 0;
        for (int l = 0; l < nLayer; l++) {
            if (!tr.v[l]) continue;
            int li = perLayer ? l : 0;
            double sum = tr.q[l][0] + tr.q[l][1] + tr.q[l][2];
            int b1 = binOf(sum, max1, nB1);
            S1 += std::log(p1[tf][0][li][b1] / p1[tf][1][li][b1]);
            for (int s = 0; s < nSlice; s++) {
                int b3 = binOf(tr.q[l][s], max3, nB3);
                S3 += std::log(p3[tf][0][li][s][b3] / p3[tf][1][li][s][b3]);
            }
            sumQ += sum;
        }
        double Sm = sumQ / tr.nv;
        double sc[3] = {S1, S3, Sm};
        for (int m = 0; m < 3; m++) (tr.sp == 0 ? sE[m] : sP[m]).push_back(sc[m]);
        hLe[0][tr.sp]->Fill(1.0 / (1.0 + std::exp(-S1)));
        hLe[1][tr.sp]->Fill(1.0 / (1.0 + std::exp(-S3)));
    }

    // ---------- 4. numbers at 90% electron efficiency ----------
    printf("\n%-14s %12s %10s\n", "method", "pion eff", "rejection");
    for (int m = 0; m < 3; m++) {
        double ep = pionEffAt(sE[m], sP[m], 0.90);
        double err = std::sqrt(ep * (1 - ep) / sP[m].size());
        printf("%-14s %6.3f+-%.3f %10.2f\n", mname[m], ep, err, ep > 0 ? 1.0 / ep : 0.0);
    }
    printf("(Felix, Run 3, 1.0-2.0 GeV/c, negative particles: LQ1D rejection 2.28, LQ3D rejection 3.36)\n");
    printf("(%zu electrons and %zu pions in the test sample; the errors cover only the pion counting)\n",
           sE[0].size(), sP[0].size());

    // ---------- 5. plots ----------
    TCanvas* c = new TCanvas("cL", "TRD likelihood (Felix definition)", 1800, 600);
    c->Divide(3, 1);

    for (int m = 0; m < 2; m++) {
        c->cd(m + 1);
        gPad->SetLogy();
        TH1D* e = hLe[m][0];
        TH1D* p = hLe[m][1];
        if (e->Integral() > 0) e->Scale(1.0 / e->Integral());
        if (p->Integral() > 0) p->Scale(1.0 / p->Integral());
        e->SetLineColor(kRed + 1);  e->SetLineWidth(2);
        p->SetLineColor(kBlue + 1); p->SetLineWidth(2);
        e->SetMaximum(3 * std::max(e->GetMaximum(), p->GetMaximum()));
        e->SetMinimum(1e-3);
        e->Draw("hist");
        p->Draw("hist same");
        TLegend* lg = new TLegend(0.40, 0.75, 0.70, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(e, Form("Electron (%zu)", sE[m].size()), "l");
        lg->AddEntry(p, Form("Pion (%zu)", sP[m].size()), "l");
        lg->Draw();
    }

    c->cd(3);
    const int col[3] = {kRed + 1, kBlue + 1, kGray + 2};
    TGraph* g[3];
    for (int m = 0; m < 3; m++) {
        g[m] = new TGraph();
        for (double ee = 0.70; ee <= 0.9801; ee += 0.01) {
            double ep = pionEffAt(sE[m], sP[m], ee);
            g[m]->SetPoint(g[m]->GetN(), ee, ep);
        }
        g[m]->SetLineColor(col[m]); g[m]->SetLineWidth(2);
    }
    g[0]->SetTitle("Pion efficiency vs electron efficiency;electron efficiency;pion efficiency");
    g[0]->GetXaxis()->SetLimits(0.69, 0.99);
    g[0]->SetMinimum(0); g[0]->SetMaximum(1);
    g[0]->Draw("AL");
    g[1]->Draw("L same");
    g[2]->Draw("L same");
    TLegend* lg3 = new TLegend(0.15, 0.68, 0.50, 0.88);
    lg3->SetBorderSize(0); lg3->SetFillStyle(0);
    for (int m = 0; m < 3; m++) lg3->AddEntry(g[m], mname[m], "l");
    lg3->Draw();

    c->SaveAs("plots/trd_likelihood_felix.png");
}