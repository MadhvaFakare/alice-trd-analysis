// trd_likelihood.C
//
// TRD electron likelihood (LQ1D and LQ3D) as in Felix Schlepper's master thesis,
// "Electron identification with the ALICE TRD in LHC Run 3 and 4" (Heidelberg 2023),
// Section 5.3, applied to the flat tree trdpid_tree.root.
//
// ---------------------------------------------------------------------------
// THE METHOD (thesis Eq. 5.3 to 5.7)
//
// 1. Reference distributions. For each momentum interval, take the charge of all
//    identified electrons and pions and make the normalised distributions P(Q|e)
//    and P(Q|pi). Felix smooths them with a fit (Landau convolved with a Gaussian,
//    his Sections 4.2.2 and 5.2.3).
//
// 2. Lookup table (LUT). Turn the two distributions into one curve:
//        L(Q) = P(Q|e) / ( P(Q|e) + P(Q|pi) )                                 (Eq. 5.3, 5.4)
//    L is near 0 where pions dominate and near 1 where electrons dominate.
//    Below a low-charge cutoff the fits are unreliable, so the LUT is fixed there.
//
// 3. LQ1D: one charge per layer, Q = Q0+Q1+Q2, one LUT.
//    LQ3D: three charges per layer, Q0, Q1 and Q2, one LUT each.
//
// 4. Combine all layers (and for LQ3D all slices) of a track:
//        L_e = prod(L) / ( prod(L) + prod(1 - L) )                            (Eq. 5.5, 5.7)
//    This is the same as L_e = 1 / (1 + exp(-S)) with S = sum of ln( L / (1-L) ).
//    The macro adds up S instead of multiplying, so nothing underflows to zero.
//
// 5. Performance: choose the L_e threshold that keeps 90% of the electrons and
//    count the fraction of pions above it (pion efficiency). Rejection = 1 / that.
//    Felix gets, for 1.0-2.0 GeV/c, negative particles: LQ1D rejection 2.28, LQ3D 3.36
//    (pp collisions at 13.6 TeV, so not directly comparable with Pb-Pb).
//
// ---------------------------------------------------------------------------
// WHERE THIS MACRO DIFFERS FROM THE THESIS, AND WHY
//
//  * Samples: Felix selected electrons and pions with TPC and TOF and says that
//    contamination was his biggest problem. His outlook proposes V0 samples
//    (gamma -> e+e-, K0s -> pi+pi-). Your fSpecies labels are V0 tags, so they are kept.
//  * Momentum: the tree has fPt, not p, and fPt has no charge sign. So the intervals
//    are in pT, and there is no split into positive and negative tracks
//    (Felix splits them because of the E x B effect).
//  * One LUT shared by all layers by default (perLayer = false). Felix makes one per
//    layer, but his Figure 5.7 shows they are very similar, and you have few electrons.
//  * Low-charge cutoff, gentle version: below the dip of the LUT (searched below the
//    electron peak), the LUT keeps the dip value instead of rising again. Felix sets the
//    LUT to 0 there. Every LUT value is also kept between 0.001 and 0.999, so a single
//    factor can never force L_e to exactly 0 or 1.
//  * Two LUT versions are always built and compared:
//      fit  : P(Q|e) and P(Q|pi) from the Landau x Gauss fits (Felix's choice)
//      hist : P(Q|e) and P(Q|pi) straight from the histograms, with a pseudo-count
//    The fit tails of both species have the same shape, so the fit LUT goes flat at high
//    charge. The histogram LUT follows the data there but is noisier.
//  * Cross validation: tracks from even DF folders are scored with LUTs built from odd
//    DF folders and the other way round. Otherwise each track would help build the LUT
//    that scores it, and the rejection would look better than it is.
//  * Errors come from a bootstrap: the electron and pion test samples are resampled
//    200 times. This includes the uncertainty of the 90% threshold. The same resamples
//    give the error on the DIFFERENCE between two methods, which is the number to look
//    at before saying one method is better.
//  * The correlation between Q0, Q1 and Q2 within a layer is printed. LQ3D multiplies
//    the three slices as if they were independent; strong correlation breaks that.
//
// A layer is valid only if Q0, Q1 and Q2 are all > 0 (-1 marks "no tracklet").
//
// Usage (from the macros folder):
//   root -l trd_likelihood.C
//   root -l 'trd_likelihood.C("/path/trdpid_tree.root", true)'        // one LUT per layer
//   root -l 'trd_likelihood.C("/path/trdpid_tree.root", false, 0)'    // layer 0 only
// ---------------------------------------------------------------------------

// ======================= settings =======================
const int    kNLayer = 6;
const int    kNSlice = 3;

// variables: 0 = Q0+Q1+Q2 (LQ1D), 1 = Q0, 2 = Q1, 3 = Q2 (LQ3D)
const int    kNVar = 4;
const char*  kVarName[kNVar] = {"Q0+Q1+Q2", "Q0", "Q1", "Q2"};
const double kVarMax[kNVar]  = {320., 128., 128., 64.}; // Q0,Q1: 7 bit (<128), Q2: 6 bit (<64)
const int    kVarBins[kNVar] = {32, 32, 32, 32};

// pT intervals, Felix's momentum intervals (he uses p)
const int    kNPt = 4;
const double kPtEdge[kNPt + 1] = {0.2, 1.0, 2.0, 3.0, 8.0};

const int    kMinTrainE     = 20;    // electrons per fold needed to use a pT interval
const int    kMinFitEntries = 50;    // entries needed to try a fit
const double kAlpha         = 0.5;   // pseudo-count per bin for the histogram pdfs
const double kLutFloor      = 1e-3;  // LUT values kept inside [floor, 1 - floor]
const bool   kHoldDip       = true;  // low-charge cutoff: below the dip, keep the dip value
const int    kNBoot         = 200;   // bootstrap samples for the errors
const double kEffE          = 0.90;  // working point: electron efficiency

// LUT kinds
const int    kNKind = 2;
const char*  kKindName[kNKind] = {"fit", "hist"};

// ======================= data structures =======================

// one tagged track
struct LqTrack {
    int    sp;                  // 0 electron, 1 pion
    int    fold;                // 0 or 1 (even or odd DF folder)
    int    ptBin;               // which pT interval
    int    nValid;              // number of valid layers
    bool   valid[kNLayer];
    double q[kNLayer][kNSlice];

    // value of variable v in layer l
    double var(int l, int v) const {
        if (v == 0) return q[l][0] + q[l][1] + q[l][2];
        return q[l][v - 1];
    }
};

// one reference distribution P(Q|species) for one variable
struct LqPdf {
    TH1D*  h     = nullptr;  // raw counts
    TH1D*  hDens = nullptr;  // density from the histogram, with pseudo-count
    TF1*   f     = nullptr;  // fit, nullptr if the fit was not tried or failed
    double norm  = 1.;       // N * bin width: turns the fitted counts into a density

    // density at x; useFit = false, or no fit available, gives the histogram density
    double eval(double x, bool useFit) const {
        double d = (useFit && f) ? f->Eval(x) / norm : hDens->Interpolate(x);
        return d > 1e-12 ? d : 1e-12;
    }
};

// ======================= helper functions =======================

// Landau convolved with a Gaussian (from the ROOT tutorial langaus.C)
// par[0] Landau width, par[1] most probable value, par[2] area, par[3] Gaussian sigma
double lqLangau(double* x, double* par)
{
    const double invsq2pi = 0.3989422804014;
    const double mpshift  = -0.22278298;
    const double np = 100.0, sc = 5.0;

    double mpc  = par[1] - mpshift * par[0];
    double xlow = x[0] - sc * par[3];
    double xupp = x[0] + sc * par[3];
    double step = (xupp - xlow) / np;
    double sum  = 0;
    for (double i = 1.0; i <= np / 2; i++) {
        double xx = xlow + (i - 0.5) * step;
        sum += TMath::Landau(xx, mpc, par[0]) / par[0] * TMath::Gaus(x[0], xx, par[3]);
        xx = xupp - (i - 0.5) * step;
        sum += TMath::Landau(xx, mpc, par[0]) / par[0] * TMath::Gaus(x[0], xx, par[3]);
    }
    return par[2] * step * sum * invsq2pi / par[3];
}

// build both densities of one pdf: histogram always, fit if possible
// returns true if the fit converged
bool lqBuildPdf(LqPdf& p)
{
    TH1D* h   = p.h;
    int    nb = h->GetNbinsX();
    double bw = h->GetBinWidth(1);
    double N  = h->Integral();

    p.hDens = (TH1D*)h->Clone(Form("%s_dens", h->GetName()));
    for (int b = 1; b <= nb; b++)
        p.hDens->SetBinContent(b, (h->GetBinContent(b) + kAlpha) / ((N + kAlpha * nb) * bw));

    p.f = nullptr;
    if (N < kMinFitEntries) return false;

    double xmax = h->GetXaxis()->GetXmax();
    double mpv  = h->GetBinCenter(h->GetMaximumBin());
    double rms  = h->GetRMS();
    TF1* f = new TF1(Form("%s_fit", h->GetName()), lqLangau, 0, xmax, 4);
    f->SetParameters(0.2 * rms, mpv, N * bw, 0.2 * rms);
    f->SetParLimits(0, 0.01, 5 * rms);
    f->SetParLimits(1, 0, xmax);
    f->SetParLimits(2, 0.2 * N * bw, 5 * N * bw);
    f->SetParLimits(3, 0.01, 5 * rms);

    // Q quiet, L Poisson likelihood (good for few entries), N do not draw, S return result
    TFitResultPtr r = h->Fit(f, "QLNS");
    if ((int)r == 0 && r->IsValid()) {
        p.f = f;
        p.norm = N * bw;
        return true;
    }
    delete f;
    return false;
}

// lookup table L(Q) = P(Q|e) / (P(Q|e) + P(Q|pi)), one bin per charge unit
TH1D* lqBuildLut(const LqPdf& e, const LqPdf& pi, bool useFit, double xmax, const char* name)
{
    int nb = (int)xmax;
    TH1D* lut = new TH1D(name, ";charge (a.u.);L_{e}", nb, 0, xmax);
    for (int b = 1; b <= nb; b++) {
        double x  = lut->GetBinCenter(b);
        double pe = e.eval(x, useFit);
        double pp = pi.eval(x, useFit);
        double L  = pe / (pe + pp);
        if (L < kLutFloor)     L = kLutFloor;
        if (L > 1 - kLutFloor) L = 1 - kLutFloor;
        lut->SetBinContent(b, L);
    }

    // low-charge cutoff (Felix 5.3.1), gentle version:
    // find the lowest LUT value below the electron peak, and keep that value for all
    // smaller charges, so the LUT cannot rise again where there is no data
    if (kHoldDip) {
        double ePeak = e.h->GetBinCenter(e.h->GetMaximumBin());
        int bMin = 1;
        for (int b = 1; b <= nb && lut->GetBinCenter(b) < ePeak; b++)
            if (lut->GetBinContent(b) < lut->GetBinContent(bMin)) bMin = b;
        for (int b = 1; b < bMin; b++) lut->SetBinContent(b, lut->GetBinContent(bMin));
    }
    return lut;
}

double lqLookup(TH1D* lut, double x)
{
    int b = lut->FindFixBin(x);
    if (b < 1) b = 1;
    if (b > lut->GetNbinsX()) b = lut->GetNbinsX();
    return lut->GetBinContent(b);
}

// pion efficiency at electron efficiency effE (higher score = more electron-like)
double lqPionEff(std::vector<double> e, const std::vector<double>& p, double effE)
{
    if (e.empty() || p.empty()) return -1;
    std::sort(e.begin(), e.end());
    size_t k = (size_t)std::floor((1.0 - effE) * e.size());
    if (k >= e.size()) k = e.size() - 1;
    double cut = e[k];                       // at least effE of the electrons are >= cut
    long n = 0;
    for (double s : p) if (s >= cut) n++;
    return (double)n / p.size();
}

// standard deviation of a list of numbers
double lqSd(const std::vector<double>& x)
{
    if (x.size() < 2) return 0;
    double s = 0, s2 = 0;
    for (double v : x) { s += v; s2 += v * v; }
    double m = s / x.size();
    return std::sqrt(std::max(0.0, s2 / x.size() - m * m));
}

// Pearson correlation from sums
double lqCorr(double n, double sx, double sy, double sxx, double syy, double sxy)
{
    double cov = sxy / n - (sx / n) * (sy / n);
    double vx  = sxx / n - (sx / n) * (sx / n);
    double vy  = syy / n - (sy / n) * (sy / n);
    return (vx > 0 && vy > 0) ? cov / std::sqrt(vx * vy) : 0;
}

// ======================= main =======================

void trd_likelihood(const char* fname = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root",
                    bool perLayer = false,  // one LUT per layer (Felix) or one shared LUT
                    int  onlyLayer = -1,    // >= 0: use only this layer
                    int  minLayers = 4)     // valid layers a track needs
{
    TFile* fin = TFile::Open(fname);
    if (!fin || fin->IsZombie()) { printf("Cannot open %s\n", fname); return; }
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.06);

    if (onlyLayer >= 0) minLayers = 1;
    const int nGroup = perLayer ? kNLayer : 1;   // how many LUTs per variable

    // ---------- 1. read the tagged tracks ----------
    std::vector<LqTrack> trk;
    int dfIndex = 0;
    long nOutsidePt = 0, nFewLayers = 0;

    TIter next(fin->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;
        TDirectory* dir = (TDirectory*)key->ReadObj();
        TTree* t = (TTree*)dir->Get("trdpid");
        if (!t) continue;

        TLeaf* lSp = t->GetLeaf("fSpecies");
        TLeaf* lPt = t->GetLeaf("fPt");
        TLeaf* lQ[kNSlice][kNLayer];
        bool ok = (lSp && lPt);
        for (int s = 0; s < kNSlice && ok; s++)
            for (int l = 0; l < kNLayer; l++) {
                lQ[s][l] = t->GetLeaf(Form("fQ%d%d", s, l));
                if (!lQ[s][l]) { ok = false; break; }
            }
        if (!ok) continue;

        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("fSpecies", 1);
        t->SetBranchStatus("fPt", 1);
        for (int s = 0; s < kNSlice; s++)
            for (int l = 0; l < kNLayer; l++)
                t->SetBranchStatus(Form("fQ%d%d", s, l), 1);

        int fold = dfIndex % 2;
        dfIndex++;

        for (Long64_t i = 0; i < t->GetEntries(); i++) {
            t->GetEntry(i);
            int sp = (int)lSp->GetValue();
            if (sp != 0 && sp != 1) continue;            // unknown (-1)

            double pt = lPt->GetValue();
            int ptBin = -1;
            for (int b = 0; b < kNPt; b++)
                if (pt >= kPtEdge[b] && pt < kPtEdge[b + 1]) ptBin = b;
            if (ptBin < 0) { nOutsidePt++; continue; }

            LqTrack tr;
            tr.sp = sp; tr.fold = fold; tr.ptBin = ptBin; tr.nValid = 0;
            for (int l = 0; l < kNLayer; l++) {
                for (int s = 0; s < kNSlice; s++) tr.q[l][s] = lQ[s][l]->GetValue();
                tr.valid[l] = tr.q[l][0] > 0 && tr.q[l][1] > 0 && tr.q[l][2] > 0
                              && (onlyLayer < 0 || l == onlyLayer);
                if (tr.valid[l]) tr.nValid++;
            }
            if (tr.nValid < minLayers) { nFewLayers++; continue; }
            trk.push_back(tr);
        }
    }

    // counts per fold, pT interval, species
    long cnt[2][kNPt][2] = {};
    for (const auto& tr : trk) cnt[tr.fold][tr.ptBin][tr.sp]++;

    printf("\n=== Tracks ===\n"
           "DF folders: %d   LUTs: %s   layers: %s   min valid layers: %d\n"
           "outside pT intervals: %ld   too few valid layers: %ld   used: %zu\n\n"
           "%-18s %18s %18s\n", dfIndex,
           perLayer ? "one per layer" : "shared by all layers",
           onlyLayer >= 0 ? Form("only layer %d", onlyLayer) : "all",
           minLayers, nOutsidePt, nFewLayers, trk.size(),
           "pT interval", "fold 0 (e / pi)", "fold 1 (e / pi)");

    bool usable[kNPt];
    int  nUsable = 0;
    for (int b = 0; b < kNPt; b++) {
        usable[b] = cnt[0][b][0] >= kMinTrainE && cnt[1][b][0] >= kMinTrainE
                 && cnt[0][b][1] >= kMinTrainE && cnt[1][b][1] >= kMinTrainE;
        if (usable[b]) nUsable++;
        printf("%4.1f - %4.1f GeV/c %10ld / %-6ld %10ld / %-6ld %s\n",
               kPtEdge[b], kPtEdge[b + 1], cnt[0][b][0], cnt[0][b][1],
               cnt[1][b][0], cnt[1][b][1], usable[b] ? "" : "<- too few, skipped");
    }
    if (nUsable == 0) { printf("No pT interval has enough tracks.\n"); return; }

    // ---------- 2. correlation of the slices within a layer ----------
    // sums per pT interval and species: index 0,1,2 = Q0,Q1,Q2
    double cN[kNPt][2] = {}, cS[kNPt][2][3] = {}, cSS[kNPt][2][3] = {}, cSP[kNPt][2][3] = {};
    const int pairA[3] = {0, 0, 1}, pairB[3] = {1, 2, 2};
    for (const auto& tr : trk)
        for (int l = 0; l < kNLayer; l++) {
            if (!tr.valid[l]) continue;
            int b = tr.ptBin, s = tr.sp;
            cN[b][s]++;
            for (int k = 0; k < 3; k++) {
                cS[b][s][k]  += tr.q[l][k];
                cSS[b][s][k] += tr.q[l][k] * tr.q[l][k];
                cSP[b][s][k] += tr.q[l][pairA[k]] * tr.q[l][pairB[k]];
            }
        }
    printf("\n=== Correlation of the slices within a layer (LQ3D assumes 0) ===\n"
           "%-18s %-9s %8s %8s %8s\n", "pT interval", "species", "Q0-Q1", "Q0-Q2", "Q1-Q2");
    for (int b = 0; b < kNPt; b++) {
        if (!usable[b]) continue;
        for (int s = 0; s < 2; s++) {
            printf("%4.1f - %4.1f GeV/c %-9s", kPtEdge[b], kPtEdge[b + 1], s == 0 ? "electron" : "pion");
            for (int k = 0; k < 3; k++) {
                int a = pairA[k], c = pairB[k];
                printf(" %8.2f", lqCorr(cN[b][s], cS[b][s][a], cS[b][s][c],
                                        cSS[b][s][a], cSS[b][s][c], cSP[b][s][k]));
            }
            printf("\n");
        }
    }

    // ---------- 3. reference distributions per fold ----------
    // pdf[fold][ptBin][group][variable][species]
    static LqPdf pdf[2][kNPt][kNLayer][kNVar][2];
    // lut[kind][fold][ptBin][group][variable], kind 0 = fit, 1 = hist
    static TH1D* lut[kNKind][2][kNPt][kNLayer][kNVar];

    for (int fo = 0; fo < 2; fo++)
        for (int b = 0; b < kNPt; b++)
            for (int g = 0; g < nGroup; g++)
                for (int v = 0; v < kNVar; v++)
                    for (int s = 0; s < 2; s++) {
                        pdf[fo][b][g][v][s] = LqPdf();
                        pdf[fo][b][g][v][s].h = new TH1D(Form("lqh_f%d_b%d_g%d_v%d_s%d", fo, b, g, v, s),
                                                         "", kVarBins[v], 0, kVarMax[v]);
                    }

    for (const auto& tr : trk)
        for (int l = 0; l < kNLayer; l++) {
            if (!tr.valid[l]) continue;
            int g = perLayer ? l : 0;
            for (int v = 0; v < kNVar; v++)
                pdf[tr.fold][tr.ptBin][g][v][tr.sp].h->Fill(tr.var(l, v));
        }

    // ---------- 4. fits and lookup tables ----------
    int nFitTry = 0, nFitOk = 0;
    for (int fo = 0; fo < 2; fo++)
        for (int b = 0; b < kNPt; b++) {
            if (!usable[b]) continue;
            for (int g = 0; g < nGroup; g++)
                for (int v = 0; v < kNVar; v++) {
                    for (int s = 0; s < 2; s++) {
                        nFitTry++;
                        if (lqBuildPdf(pdf[fo][b][g][v][s])) nFitOk++;
                    }
                    for (int k = 0; k < kNKind; k++)
                        lut[k][fo][b][g][v] = lqBuildLut(pdf[fo][b][g][v][0], pdf[fo][b][g][v][1], k == 0,
                                                         kVarMax[v],
                                                         Form("lqlut_k%d_f%d_b%d_g%d_v%d", k, fo, b, g, v));
                }
        }
    printf("\nFits: %d of %d converged (a failed fit falls back to the histogram)\n", nFitOk, nFitTry);

    // ---------- 5. score every track with the LUTs of the other fold ----------
    const int   nMeth = 5;
    const char* methName[nMeth] = {"LQ1D fit", "LQ3D fit", "LQ1D hist", "LQ3D hist", "mean charge"};
    std::vector<double> score[kNPt][nMeth][2];   // same track order in every method
    TH1D* hLe[kNPt][2][2];                       // L_e of the fit methods, [ptBin][LQ1D/LQ3D][species]
    for (int b = 0; b < kNPt; b++)
        for (int m = 0; m < 2; m++)
            for (int s = 0; s < 2; s++)
                hLe[b][m][s] = new TH1D(Form("lqLe_b%d_m%d_s%d", b, m, s),
                                        Form("%s, %.1f-%.1f GeV/c;L_{e};fraction of tracks",
                                             methName[m], kPtEdge[b], kPtEdge[b + 1]),
                                        50, 0, 1);

    for (const auto& tr : trk) {
        int b = tr.ptBin;
        if (!usable[b]) continue;
        int fo = 1 - tr.fold;                    // LUTs from the other fold
        double S1[kNKind] = {0, 0}, S3[kNKind] = {0, 0}, sumQ = 0;
        for (int l = 0; l < kNLayer; l++) {
            if (!tr.valid[l]) continue;
            int g = perLayer ? l : 0;
            for (int k = 0; k < kNKind; k++) {
                double L = lqLookup(lut[k][fo][b][g][0], tr.var(l, 0));
                S1[k] += std::log(L / (1 - L));
                for (int v = 1; v < kNVar; v++) {
                    L = lqLookup(lut[k][fo][b][g][v], tr.var(l, v));
                    S3[k] += std::log(L / (1 - L));
                }
            }
            sumQ += tr.var(l, 0);
        }
        double sc[nMeth] = {S1[0], S3[0], S1[1], S3[1], sumQ / tr.nValid};
        for (int m = 0; m < nMeth; m++) score[b][m][tr.sp].push_back(sc[m]);
        hLe[b][0][tr.sp]->Fill(1.0 / (1.0 + std::exp(-S1[0])));
        hLe[b][1][tr.sp]->Fill(1.0 / (1.0 + std::exp(-S3[0])));
    }

    // ---------- 6. results with bootstrap errors ----------
    printf("\n=== Pion efficiency at %.0f%% electron efficiency (bootstrap errors) ===\n", 100 * kEffE);
    for (int b = 0; b < kNPt; b++) {
        if (!usable[b]) continue;
        size_t nE = score[b][0][0].size(), nP = score[b][0][1].size();

        // resample the same tracks for every method, so differences get correct errors
        TRandom3 rnd(12345 + b);
        std::vector<double> boot[nMeth];
        std::vector<size_t> ie(nE), ip(nP);
        std::vector<double> eb(nE), pb(nP);
        for (int k = 0; k < kNBoot; k++) {
            for (size_t i = 0; i < nE; i++) ie[i] = rnd.Integer(nE);
            for (size_t i = 0; i < nP; i++) ip[i] = rnd.Integer(nP);
            for (int m = 0; m < nMeth; m++) {
                for (size_t i = 0; i < nE; i++) eb[i] = score[b][m][0][ie[i]];
                for (size_t i = 0; i < nP; i++) pb[i] = score[b][m][1][ip[i]];
                boot[m].push_back(lqPionEff(eb, pb, kEffE));
            }
        }

        printf("\npT %.1f - %.1f GeV/c: %zu electrons, %zu pions\n",
               kPtEdge[b], kPtEdge[b + 1], nE, nP);
        printf("  %-12s %18s %10s\n", "method", "pion efficiency", "rejection");
        double val[nMeth];
        for (int m = 0; m < nMeth; m++) {
            val[m] = lqPionEff(score[b][m][0], score[b][m][1], kEffE);
            printf("  %-12s %9.3f +- %.3f %10.2f\n", methName[m], val[m], lqSd(boot[m]),
                   val[m] > 0 ? 1.0 / val[m] : 0.0);
        }

        // differences: significant only if well outside the error
        const int dA[4] = {1, 3, 2, 3}, dB[4] = {0, 2, 0, 1};
        printf("  differences (negative = first method better):\n");
        for (int d = 0; d < 4; d++) {
            std::vector<double> diff(kNBoot);
            for (int k = 0; k < kNBoot; k++) diff[k] = boot[dA[d]][k] - boot[dB[d]][k];
            double dv = val[dA[d]] - val[dB[d]], de = lqSd(diff);
            printf("    %-10s - %-10s %+7.3f +- %.3f  (%.1f sigma)\n",
                   methName[dA[d]], methName[dB[d]], dv, de, de > 0 ? std::fabs(dv) / de : 0.0);
        }
    }
    printf("\nFelix (pp 13.6 TeV, TPC+TOF samples, p not pT, negative tracks, 1.0-2.0 GeV/c): "
           "LQ1D rejection 2.28, LQ3D 3.36\n");

    // ---------- 7. plots ----------
    gSystem->mkdir("plots", true);
    const int colE = kRed + 1, colP = kBlue + 1;
    const int layerCol[kNLayer] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1};

    // 7a. reference distributions with fits (fold 0, layer group 0)
    TCanvas* cPdf = new TCanvas("lqPdf", "Reference distributions (fold 0)", 1800, 450 * nUsable);
    cPdf->Divide(kNVar, nUsable);
    int row = 0;
    for (int b = 0; b < kNPt; b++) {
        if (!usable[b]) continue;
        for (int v = 0; v < kNVar; v++) {
            cPdf->cd(row * kNVar + v + 1);
            double ymax = 0;
            TH1D* hd[2];
            TGraph* gf[2];
            for (int s = 0; s < 2; s++) {
                const LqPdf& p = pdf[0][b][0][v][s];
                hd[s] = (TH1D*)p.h->Clone(Form("lqd_b%d_v%d_s%d", b, v, s));
                double N = hd[s]->Integral();
                if (N > 0) hd[s]->Scale(1.0 / (N * hd[s]->GetBinWidth(1)));
                hd[s]->SetLineColor(s == 0 ? colE : colP);
                hd[s]->SetLineWidth(2);
                ymax = std::max(ymax, hd[s]->GetMaximum());
                gf[s] = nullptr;
                if (p.f) {
                    gf[s] = new TGraph();
                    for (int i = 0; i <= 200; i++) {
                        double x = kVarMax[v] * i / 200.;
                        gf[s]->SetPoint(i, x, p.f->Eval(x) / p.norm);
                    }
                    gf[s]->SetLineColor(s == 0 ? colE : colP);
                    gf[s]->SetLineStyle(2);
                    gf[s]->SetLineWidth(2);
                }
            }
            hd[0]->SetTitle(Form("%s, %.1f-%.1f GeV/c;charge (a.u.);density",
                                 kVarName[v], kPtEdge[b], kPtEdge[b + 1]));
            hd[0]->SetMaximum(1.2 * ymax);
            hd[0]->Draw("hist");
            hd[1]->Draw("hist same");
            for (int s = 0; s < 2; s++) if (gf[s]) gf[s]->Draw("L same");
            if (row == 0 && v == 0) {
                TLegend* lg = new TLegend(0.55, 0.70, 0.88, 0.88);
                lg->SetBorderSize(0); lg->SetFillStyle(0);
                lg->AddEntry(hd[0], "Electron", "l");
                lg->AddEntry(hd[1], "Pion", "l");
                lg->Draw();
            }
        }
        row++;
    }
    cPdf->SaveAs("plots/trd_likelihood_pdfs.png");

    // 7b. lookup tables (fold 0): fit solid, histogram dashed. Like Felix's Figures 5.7 and 5.9
    TCanvas* cLut = new TCanvas("lqLut", "Lookup tables (fold 0)", 1800, 450 * nUsable);
    cLut->Divide(kNVar, nUsable);
    row = 0;
    for (int b = 0; b < kNPt; b++) {
        if (!usable[b]) continue;
        for (int v = 0; v < kNVar; v++) {
            cLut->cd(row * kNVar + v + 1);
            bool first = true;
            for (int g = 0; g < nGroup; g++)
                for (int k = 0; k < kNKind; k++) {
                    TH1D* h = lut[k][0][b][g][v];
                    h->SetTitle(Form("LUT %s, %.1f-%.1f GeV/c;charge (a.u.);L_{e}",
                                     kVarName[v], kPtEdge[b], kPtEdge[b + 1]));
                    h->SetLineColor(layerCol[g]);
                    h->SetLineStyle(k == 0 ? 1 : 2);
                    h->SetLineWidth(2);
                    h->SetMinimum(0);
                    h->SetMaximum(1.05);
                    h->Draw(first ? "L" : "L same");
                    first = false;
                }
            if (row == 0 && v == 0) {
                TLegend* lg = new TLegend(0.55, 0.15, 0.88, 0.40);
                lg->SetBorderSize(0); lg->SetFillStyle(0);
                lg->AddEntry(lut[0][0][b][0][v], "fit", "l");
                lg->AddEntry(lut[1][0][b][0][v], "histogram", "l");
                if (perLayer)
                    for (int g = 1; g < nGroup; g++)
                        lg->AddEntry(lut[0][0][b][g][v], Form("layer %d", g), "l");
                lg->Draw();
            }
        }
        row++;
    }
    cLut->SaveAs("plots/trd_likelihood_luts.png");

    // 7c. L_e distributions (fit LUTs) and pion vs electron efficiency for all methods
    TCanvas* cRes = new TCanvas("lqRes", "Electron likelihood", 1800, 500 * nUsable);
    cRes->Divide(3, nUsable);
    const int methCol[nMeth]   = {kRed + 1, kBlue + 1, kRed + 1, kBlue + 1, kGray + 2};
    const int methStyle[nMeth] = {1, 1, 2, 2, 1};
    row = 0;
    for (int b = 0; b < kNPt; b++) {
        if (!usable[b]) continue;
        for (int m = 0; m < 2; m++) {
            cRes->cd(row * 3 + m + 1);
            gPad->SetLogy();
            TH1D* e = hLe[b][m][0];
            TH1D* p = hLe[b][m][1];
            if (e->Integral() > 0) e->Scale(1.0 / e->Integral());
            if (p->Integral() > 0) p->Scale(1.0 / p->Integral());
            e->SetLineColor(colE); e->SetLineWidth(2);
            p->SetLineColor(colP); p->SetLineWidth(2);
            e->SetMaximum(3 * std::max(e->GetMaximum(), p->GetMaximum()));
            e->SetMinimum(1e-3);
            e->Draw("hist");
            p->Draw("hist same");
            TLegend* lg = new TLegend(0.35, 0.75, 0.70, 0.88);
            lg->SetBorderSize(0); lg->SetFillStyle(0);
            lg->AddEntry(e, Form("Electron (%zu)", score[b][m][0].size()), "l");
            lg->AddEntry(p, Form("Pion (%zu)", score[b][m][1].size()), "l");
            lg->Draw();
        }
        cRes->cd(row * 3 + 3);
        TGraph* gr[nMeth];
        for (int m = 0; m < nMeth; m++) {
            gr[m] = new TGraph();
            for (int i = 0; i <= 28; i++) {
                double ee = 0.70 + 0.01 * i;
                gr[m]->SetPoint(i, ee, lqPionEff(score[b][m][0], score[b][m][1], ee));
            }
            gr[m]->SetLineColor(methCol[m]);
            gr[m]->SetLineStyle(methStyle[m]);
            gr[m]->SetLineWidth(2);
        }
        gr[0]->SetTitle(Form("%.1f-%.1f GeV/c;electron efficiency;pion efficiency",
                             kPtEdge[b], kPtEdge[b + 1]));
        gr[0]->GetXaxis()->SetLimits(0.69, 0.99);
        gr[0]->SetMinimum(0);
        gr[0]->SetMaximum(1);
        gr[0]->Draw("AL");
        for (int m = 1; m < nMeth; m++) gr[m]->Draw("L same");
        TLegend* lg = new TLegend(0.15, 0.60, 0.50, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        for (int m = 0; m < nMeth; m++) lg->AddEntry(gr[m], methName[m], "l");
        lg->Draw();
        row++;
    }
    cRes->SaveAs("plots/trd_likelihood_results.png");
}