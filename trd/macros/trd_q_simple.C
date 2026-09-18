void trd_q_simple() {

    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/trd/AnalysisResults.root");
    if (!f || f->IsZombie()) {
        printf("Error opening file\n");
        return;
    }

    // ---- get all six histograms ----
    TH1F *hQ0El = (TH1F*)f->Get("t-r-d-p-i-d-study/Q0_El");
    TH1F *hQ0Pi = (TH1F*)f->Get("t-r-d-p-i-d-study/Q0_Pi");
    TH1F *hQ1El = (TH1F*)f->Get("t-r-d-p-i-d-study/Q1_El");
    TH1F *hQ1Pi = (TH1F*)f->Get("t-r-d-p-i-d-study/Q1_Pi");
    TH1F *hQ2El = (TH1F*)f->Get("t-r-d-p-i-d-study/Q2_El");
    TH1F *hQ2Pi = (TH1F*)f->Get("t-r-d-p-i-d-study/Q2_Pi");

    if (!hQ0El || !hQ0Pi || !hQ1El || !hQ1Pi || !hQ2El || !hQ2Pi) {
        printf("Could not find one of the histograms\n");
        return;
    }

    // ---- means, read before scaling ----
    double mQ0El = hQ0El->GetMean();
    double mQ0Pi = hQ0Pi->GetMean();
    double mQ1El = hQ1El->GetMean();
    double mQ1Pi = hQ1Pi->GetMean();
    double mQ2El = hQ2El->GetMean();
    double mQ2Pi = hQ2Pi->GetMean();

    printf("Q0: electrons %.2f, pions %.2f\n", mQ0El, mQ0Pi);
    printf("Q1: electrons %.2f, pions %.2f\n", mQ1El, mQ1Pi);
    printf("Q2: electrons %.2f, pions %.2f\n", mQ2El, mQ2Pi);

    // ---- normalise so the shapes can be compared ----
    hQ0El->Scale(1.0/hQ0El->Integral());
    hQ0Pi->Scale(1.0/hQ0Pi->Integral());
    hQ1El->Scale(1.0/hQ1El->Integral());
    hQ1Pi->Scale(1.0/hQ1Pi->Integral());
    hQ2El->Scale(1.0/hQ2El->Integral());
    hQ2Pi->Scale(1.0/hQ2Pi->Integral());

    // ---- colours ----
    hQ0El->SetLineColor(kRed);   hQ0El->SetLineWidth(2);
    hQ1El->SetLineColor(kRed);   hQ1El->SetLineWidth(2);
    hQ2El->SetLineColor(kRed);   hQ2El->SetLineWidth(2);
    hQ0Pi->SetLineColor(kBlue);  hQ0Pi->SetLineWidth(2);
    hQ1Pi->SetLineColor(kBlue);  hQ1Pi->SetLineWidth(2);
    hQ2Pi->SetLineColor(kBlue);  hQ2Pi->SetLineWidth(2);

    hQ0El->SetTitle("Q0: electrons (red) vs pions (blue);Q0 (a.u.);Normalised entries");
    hQ1El->SetTitle("Q1: electrons (red) vs pions (blue);Q1 (a.u.);Normalised entries");
    hQ2El->SetTitle("Q2: electrons (red) vs pions (blue);Q2 (a.u.);Normalised entries");

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "TRD charge slices", 1500, 500);
    c->Divide(3, 1);

    c->cd(1);
    hQ0El->Draw("HIST");
    hQ0Pi->Draw("HIST SAME");

    c->cd(2);
    hQ1El->Draw("HIST");
    hQ1Pi->Draw("HIST SAME");

    c->cd(3);
    hQ2El->Draw("HIST");
    hQ2Pi->Draw("HIST SAME");

    c->SaveAs("/home/madhvafakare/alice/my-analysis/plots/trd_q_simple.png");
}