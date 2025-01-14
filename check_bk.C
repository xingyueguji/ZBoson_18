void check_bk(int opt = 1)
{
    // 1 = raw, 2 = eff, 3 = eta, 4 = eta + eff
    Int_t cenlowlimit[11] = {0, 10, 20, 30, 30, 0, 15, 50, 0, 14, 0};
    Int_t cenhighlimit[11] = {10, 20, 30, 100, 50, 15, 100, 100, 14, 100, 100};

    TString path = "./rootfile/normalized/";
    TString file;
    TString marker;
    if (opt == 1)
    {
        file = "rawfile.root";
        marker = "raw";
    }
    if (opt == 2)
    {
        file = "efffile.root";
        marker = "eff";
    }
    if (opt == 3)
    {
        file = "etacut_file.root";
        marker = "eta";
    }
    if (opt == 4)
    {
        file = "etacut_eff_file.root";
        marker = "eta+eff";
    }

    TString total_file_path = path + file;

    TFile *f1 = new TFile(total_file_path, "READ");
    TH1D *bk[5];

    for (int i = 0; i < 11; i++)
    {
        if (i < 4)
        {
            bk[i] = (TH1D *)f1->Get(Form("Normalized_mc_bk_%i", i));
        }
        if (i == 10)
        {
            bk[4] = (TH1D *)f1->Get("Normalized_mc_bk_10");
        }
    }

    for (int i = 0; i < 5; i++)
    {
        TCanvas *c1 = new TCanvas("", "", 800, 600);
        c1->cd();
        // bk[i]->Rebin(4);
        bk[i]->SetMarkerStyle(20);   // Set marker style (e.g., solid circles)
        bk[i]->SetMarkerSize(1.0);   // Set marker size
        bk[i]->SetMarkerColor(kRed); // Set marker color (use ROOT color constants)
        bk[i]->SetLineColor(kBlue);  // Set line color
        bk[i]->SetLineWidth(2);      // Set line width
        bk[i]->SetFillColor(kGreen); // Set fill color for histogram
        bk[i]->SetFillStyle(3004);   // Set fill style (e.g., hatched lines)

        bk[i]->Draw("E1"); // Draw with error bars ("E1" for errors with lines)
        c1->SaveAs(Form("./background/bk_%s_%i.png", marker.Data(), i));
        delete c1;
    }
}