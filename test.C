void test()
{
    TFile *f1 = new TFile("./rootfile/template_PbPb.root");
    TFile *f2 = new TFile("./rootfile/mc_signal.root");

    TH1D *raw = (TH1D *)f2->Get("FA_nominal_10");
    TH1D *notraw = (TH1D *)f1->Get("template_FA_nominal_21_41_10");

    cout << "nentriesraw is " << raw->GetEntries() << endl;
    cout << "nentriesnotraw is " << notraw->GetEntries() << endl;

    cout << "diff is " << notraw->GetEntries() - raw->GetEntries() << endl;

    TCanvas *c1 = new TCanvas("", "", 800, 800);
    c1->SetLogy();
    c1->cd();

    raw->Scale(1/ raw->Integral("width"));
    notraw->Scale(1/ notraw->Integral("width"));

    raw->Rebin(4);
    raw->SetMarkerStyle(22);
    raw->SetMarkerColor(kRed);
    notraw->Rebin(4);
    notraw->SetMarkerStyle(22);
    notraw->SetMarkerColor(kGreen);

    raw->Draw("P");
    notraw->Draw("P SAME");

    //raw->Draw("TEXT SAME");
    //notraw->Draw("TEXT SAME");




}