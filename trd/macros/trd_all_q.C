/*
 * Electron / pion comparison for every charge column in the tree, plus the
 * diagnostics that decide whether the numbers can be trusted.
 *
 * Slice definitions, from Table 4.1 of F. Schlepper's master's thesis
 * (each time bin is 100 ns):
 *
 *   Q0  bins  1-6   amplification peak   7 bit, saturates at 127
 *   Q1  bins 14-20  TR peak              7 bit, saturates at 127
 *   Q2  bins  7-13  plateau region       6 bit, saturates at  63
 *
 * The slices are not in time order: Q0 is earliest, Q2 is in the middle, and
 * Q1 is latest. Q1 covers the TR peak, so it should separate best.
 *
 * Charge correction, Schlepper Eq. 5.1, implemented by Kangkan as
 *
 *   Q_corr = Q / (C * sqrt(1 + tan^2(phi) + tan^2(lambda)))
 *
 * where C is the local gain factor. The square-root term normalises to the
 * track length.
 *
 * This macro reports two hardware effects, both from Table 4.1:
 *
 * Saturation:
 *   Each slice is divided by a scale factor of 16 on the MCM, and values above
 *   the range are truncated to the maximum. The high-charge tail is where the
 *   electron excess lives. If many tracks sit near the ceiling, the measured
 *   separation is systematically low.
 *
 * Low-end artefact:
 *   Schlepper notes entries at low charge in Q0 and Q1 caused by binning with
 *   that scale factor, and excludes them from his fits. PIDStudy.cxx cuts below
 *   6, which may or may not be the correct threshold.
 *
 * Usage:
 *   root -l 'macros/trd_all_q_refactored.C'
 *   root -l -b -q 'macros/trd_all_q_refactored.C' > all_q.txt
 */

/*
 * Return the pion efficiency at a target electron efficiency when cutting
 * from above. Integrating through nb + 1 includes the overflow bin.
 */
double pionEffAt(TH1F *hE, TH1F *hP, double target, double &cut)
{
    const int nb = hE->GetNbinsX();
    const double totE = hE->Integral(1, nb + 1);
    const double totP = hP->Integral(1, nb + 1);

    if (totE <= 0 || totP <= 0) {
        return -1;
    }

    for (int bin = nb + 1; bin >= 1; --bin) {
        if (hE->Integral(bin, nb + 1) / totE >= target) {
            cut = hE->GetXaxis()->GetBinLowEdge(bin);
            return hP->Integral(bin, nb + 1) / totP;
        }
    }

    return -1;
}

void trd_all_q(
    const char *inFile = "/home/madhvafakare/alice/PbPb2023pass5/trd/trdpid_tree.root")
{
    /* Chain every dataframe into one tree. */
    TFile *file = TFile::Open(inFile);
    if (!file || file->IsZombie()) {
        printf("Cannot open %s\n", inFile);
        return;
    }

    TChain *tree = new TChain("trdpid");
    TIter next(file->GetListOfKeys());
    TKey *key = nullptr;

    while ((key = static_cast<TKey *>(next()))) {
        const TString name = key->GetName();
        if (name.BeginsWith("DF_")) {
            tree->Add(TString(inFile) + "/" + name + "/trdpid");
        }
    }

    printf("Chained %lld tracks\n", tree->GetEntries());

    /*
     * Part 1: saturation check
     *
     * The quality window in PIDStudy.cxx keeps 6 <= Q < 127 for Q0 and Q1,
     * and 6 <= Q < 62 for Q2. Anything truncated by the hardware sits at the
     * ceiling and is discarded. Count how close the surviving tracks are to
     * that ceiling.
     */
    printf("\n=== PART 1: how close to the hardware ceiling? ===\n"
           "Q0 and Q1 saturate at 127 (7 bit), Q2 at 63 (6 bit).\n"
           "Corrected charges can exceed those limits because the correction divides\n"
           "by the gain and the path length.\n\n");

    printf("%-6s %-6s %10s %12s %12s %10s\n",
           "slice", "layer", "N usable", "N above 110", "N above 120", "max");

    for (int slice = 0; slice < 3; ++slice) {
        for (int layer = 0; layer < 6; ++layer) {
            const TString variable = Form("fQ%d%d", slice, layer);
            const double lowerCeilingThreshold = (slice == 2) ? 50 : 110;
            const double upperCeilingThreshold = (slice == 2) ? 58 : 120;
            const Long64_t nUsable =
                tree->GetEntries(Form("%s>0", variable.Data()));
            const Long64_t nAboveLower = tree->GetEntries(
                Form("%s>%f", variable.Data(), lowerCeilingThreshold));
            const Long64_t nAboveUpper = tree->GetEntries(
                Form("%s>%f", variable.Data(), upperCeilingThreshold));

            printf("Q%-5d %-6d %10lld %12lld %12lld %10.1f\n",
                   slice,
                   layer,
                   nUsable,
                   nAboveLower,
                   nAboveUpper,
                   tree->GetMaximum(variable));
        }
    }

    printf("\nFor Q2 the two thresholds are 50 and 58 because its ceiling is lower.\n");

    /*
     * Part 2: low-charge artefact
     *
     * Schlepper excludes low-charge entries in Q0 and Q1 caused by the
     * scale-factor binning. PIDStudy.cxx cuts below 6. Measure how much lies
     * just above that cut, where the artefact would appear.
     */
    printf("\n=== PART 2: how much sits just above the lower cut of 6? ===\n");
    printf("%-6s %10s %12s %12s\n", "slice", "N usable", "N below 10", "fraction");

    for (int slice = 0; slice < 3; ++slice) {
        const TString variable = Form("fQ%d0", slice);
        const Long64_t nUsable =
            tree->GetEntries(Form("%s>0", variable.Data()));
        const Long64_t nBelowTen = tree->GetEntries(
            Form("%s>0 && %s<10", variable.Data(), variable.Data()));
        const double lowFraction =
            nUsable > 0 ? 100.0 * nBelowTen / nUsable : 0.0;

        printf("Q%-5d %10lld %12lld %11.1f%%\n",
               slice, nUsable, nBelowTen, lowFraction);
    }

    printf("(layer 0 only. Schlepper excluded this region from his fits.)\n");

    /* Part 3: separation table. */
    constexpr double xMax = 600.0;
    constexpr int nBins = 300;

    gStyle->SetOptStat(0);
    TCanvas *canvas =
        new TCanvas("cAllQ", "All charge columns", 1800, 950);
    canvas->Divide(6, 3);

    const char *description[3] = {"amplification", "TR peak", "plateau"};

    printf("\n=== PART 3: all 18 charge columns, electrons vs pions ===\n"
           "(charge > 0 throughout, so the -1 no-tracklet markers are excluded)\n\n");
    printf("%-5s %-14s %-6s %9s %6s %9s %7s %10s %8s %11s %6s\n",
           "slice",
           "covers",
           "layer",
           "mean e",
           "N e",
           "mean pi",
           "N pi",
           "sep.power",
           "Q cut",
           "eff_pi@90%",
           "over");

    double bestSeparation = -1;
    TString bestVariable;

    for (int slice = 0; slice < 3; ++slice) {
        for (int layer = 0; layer < 6; ++layer) {
            const TString variable = Form("fQ%d%d", slice, layer);
            const TString title = Form(
                "Q%d layer %d;charge (a.u.);normalised", slice, layer);

            TH1F *electronHistogram = new TH1F(
                Form("hE_%d_%d", slice, layer), title, nBins, 0, xMax);
            TH1F *pionHistogram = new TH1F(
                Form("hP_%d_%d", slice, layer), title, nBins, 0, xMax);

            /*
             * Include out-of-range entries in the mean and RMS so that
             * GetEntries() and GetMean() describe the same sample.
             */
            electronHistogram->StatOverflows(kTRUE);
            pionHistogram->StatOverflows(kTRUE);

            tree->Project(electronHistogram->GetName(),
                          variable,
                          Form("fSpecies==0 && %s>0", variable.Data()));
            tree->Project(pionHistogram->GetName(),
                          variable,
                          Form("fSpecies==1 && %s>0", variable.Data()));

            const double nElectrons = electronHistogram->GetEntries();
            const double nPions = pionHistogram->GetEntries();
            const double overflow =
                electronHistogram->GetBinContent(nBins + 1) +
                pionHistogram->GetBinContent(nBins + 1);

            if (nElectrons < 10 || nPions < 10) {
                printf(
                    "Q%-4d %-14s %-6d %9s %6.0f %9s %7.0f %10s %8s %11s %6.0f\n",
                    slice,
                    description[slice],
                    layer,
                    "-",
                    nElectrons,
                    "-",
                    nPions,
                    "-",
                    "-",
                    "too few",
                    overflow);
            } else {
                const double electronMean = electronHistogram->GetMean();
                const double pionMean = pionHistogram->GetMean();
                const double electronRms = electronHistogram->GetRMS();
                const double pionRms = pionHistogram->GetRMS();
                const double averageRms = 0.5 * (electronRms + pionRms);
                const double separation = averageRms > 0
                    ? (electronMean - pionMean) / averageRms
                    : 0;

                double chargeCut = 0;
                const double pionEfficiency = pionEffAt(
                    electronHistogram, pionHistogram, 0.90, chargeCut);

                printf(
                    "Q%-4d %-14s %-6d %9.2f %6.0f %9.2f %7.0f %10.3f %8.1f %11.3f %6.0f\n",
                    slice,
                    description[slice],
                    layer,
                    electronMean,
                    nElectrons,
                    pionMean,
                    nPions,
                    separation,
                    chargeCut,
                    pionEfficiency,
                    overflow);

                if (separation > bestSeparation) {
                    bestSeparation = separation;
                    bestVariable = variable;
                }
            }

            /* Draw the normalised electron and pion distributions. */
            canvas->cd(slice * 6 + layer + 1);
            gPad->SetLogy();

            TH1F *electrons = static_cast<TH1F *>(
                electronHistogram->Clone(Form("e_%d_%d", slice, layer)));
            TH1F *pions = static_cast<TH1F *>(
                pionHistogram->Clone(Form("p_%d_%d", slice, layer)));

            if (electrons->Integral() > 0) {
                electrons->Scale(1.0 / electrons->Integral());
            }
            if (pions->Integral() > 0) {
                pions->Scale(1.0 / pions->Integral());
            }

            electrons->SetLineColor(kBlue);
            electrons->SetLineWidth(2);
            pions->SetLineColor(kRed);
            pions->SetLineWidth(2);
            electrons->SetTitle(title);
            pions->SetTitle(title);  // Either histogram may be drawn first.

            electrons->GetXaxis()->SetRangeUser(0, 160);
            pions->GetXaxis()->SetRangeUser(0, 160);

            if (electrons->GetMaximum() > pions->GetMaximum()) {
                electrons->Draw("HIST");
                pions->Draw("HIST SAME");
            } else {
                pions->Draw("HIST");
                electrons->Draw("HIST SAME");
            }

            if (slice == 0 && layer == 0) {
                TLegend *legend = new TLegend(0.45, 0.70, 0.88, 0.88);
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                legend->AddEntry(electrons, "electrons", "l");
                legend->AddEntry(pions, "pions", "l");
                legend->Draw();
            }
        }

        printf("\n");
    }

    /* Part 4: per-track sum, which requires the tree. */
    printf("=== PART 4: per-track sum Q0+Q1+Q2, layer 0 ===\n");

    TH1F electronSum("hSE", "", nBins, 0, xMax);
    TH1F pionSum("hSP", "", nBins, 0, xMax);
    electronSum.StatOverflows(kTRUE);
    pionSum.StatOverflows(kTRUE);

    tree->Project("hSE",
                  "fQ00+fQ10+fQ20",
                  "fSpecies==0 && fQ00>0 && fQ10>0 && fQ20>0");
    tree->Project("hSP",
                  "fQ00+fQ10+fQ20",
                  "fSpecies==1 && fQ00>0 && fQ10>0 && fQ20>0");

    const double electronMean = electronSum.GetMean();
    const double pionMean = pionSum.GetMean();
    const double electronMeanError = electronSum.GetRMS() /
        TMath::Sqrt(TMath::Max(electronSum.GetEntries(), 1.0));
    const double pionMeanError = pionSum.GetRMS() /
        TMath::Sqrt(TMath::Max(pionSum.GetEntries(), 1.0));

    printf(
        "electrons: %.2f +/- %.2f  (N %.0f)\n"
        "pions    : %.2f +/- %.2f  (N %.0f)\n"
        "excess   : %.1f%%\n",
        electronMean,
        electronMeanError,
        electronSum.GetEntries(),
        pionMean,
        pionMeanError,
        pionSum.GetEntries(),
        100.0 * (electronMean - pionMean) / pionMean);

    if (bestSeparation > 0) {
        printf("\nBest separation power: %s at %.3f\n",
               bestVariable.Data(), bestSeparation);
        printf("Q1 is expected to lead because it covers the TR peak.\n");
    }

    printf("\nHow to read this:\n"
           "  sep.power  difference of means over the average RMS. Below 1 means\n"
           "             the distributions overlap heavily.\n"
           "  eff_pi@90%% fraction of pions surviving a cut that keeps 90%% of\n"
           "             electrons. Lower is better. The published ALICE design\n"
           "             goal is below 0.01 with all six chambers combined\n"
           "             (NIM A 881 (2018) 88), so one chamber will be far above.\n"
           "  over       entries beyond the histogram range, included in the mean\n"
           "             via StatOverflows. This should be zero.\n"
           "\nCaveats:\n"
           "  About 220 tagged electrons in total, so each layer holds only a few\n"
           "  dozen. Single-layer numbers are indicative, not conclusive.\n"
           "  Layer 0 is the anomalous layer in Schlepper's study, attributed to\n"
           "  secondaries raising multiplicity. Cross-check in another layer.\n"
           "  Hardware saturation removes part of the high-charge tail, which is\n"
           "  where the electron excess lives. The separation measured here is\n"
           "  therefore a lower bound. See Part 1.\n\n");

    canvas->SaveAs(
        "/home/madhvafakare/alice/my-analysis/trd/plots/trd_all_q.png");
    printf("Written: trd/plots/trd_all_q.png\n");
}
