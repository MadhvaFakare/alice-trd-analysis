void trd_q0_comparison() {
    TFile *f = TFile::Open("/home/madhvafakare/alice/PbPb2023pass5/trd/AnalysisResults.root");
    if (!f || f->IsZombie()) {
        printf("Error Opening file\n");
        return;
    }
    TH1F *hEl = (TH1F*)f->Get("t-r-d-p-i-d-study/Q0_El");
    TH1F *hPi = (TH1F*)f->Get("t-r-d-p-i-d-study/Q0_Pi");
    if (!hEl || !hPi){
        printf("Could not find the histogram\n");
        return;
    }
    hEl->Scale(1.0/hEl->Integral());
    hPi->Scale(1.0/hPi->Integral());
    TCanvas *c = new TCanvas("cQ0Comparison","TRD Q0: electron vs pions", 900,700);
    hEl->SetLineColor(kRed);
    hEl->SetLineWidth(2);
    hEl->SetTitle("TRD Q0 charge: electron vs pion (V0-tagged);Q0 (a.u);Normalized entries");

    hPi->SetLineColor(kBlue);
    hPi->SetLineWidth(2);

    gStyle->SetOptStat(0);
    hEl->Draw("HIST");
    hPi->Draw("HISt SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    TLegendEntry *e1 = leg->AddEntry(hEl, Form("Electrons (mean=%.2f, N=%.0f)",hEl->GetMean(), hEl->GetEntries()), "1");
    e1->SetTextColor(kRed);
    TLegendEntry *e2 = leg->AddEntry(hPi, Form("Pions (mean=%.2f, N=%.0f)", hPi->GetMean(), hPi->GetEntries()), "1");
    e2->SetTextColor(kBlue);
    leg->Draw();
    //c->SetLogy();
    c->SaveAs("/home/madhvafakare/alice/my-analysis/plots/trd_q0_comparison.png");
}