void A_ratio(){
	TFile *file = new TFile("./mc_bk_plot/W_histogram.root","READ");
	file->cd();
	TH1D* h_1;
	TH1D* h_2;
	TH1D* h_3;

	h_1 = (TH1D*) file->Get("data_A_24");
	h_2 = (TH1D*) file->Get("data_A_100");
	h_3 = (TH1D*) h_1->Clone();

	double data_A_normalization_factor_24 = h_1->Integral("width");
	double data_A_normalization_factor_100 = h_2->Integral("width");

	h_1->Scale(1/data_A_normalization_factor_24);
	h_2->Scale(1/data_A_normalization_factor_100);
	h_3->Scale(1/data_A_normalization_factor_24);
	cout << "Area is " << h_1->Integral("width") << "Area is " << h_2->Integral("width") << endl;


	int binnumber = h_1->GetNbinsX();
	//cout << binnumber << endl;

	/*for (int i = 1; i <= binnumber; i++){
		double ratio = 0;
		if (h_2->GetBinContent(i) != 0)ratio = h_1->GetBinContent(i) / h_2->GetBinContent(i);
		h_3->SetBinContent(i,ratio);
	}*/
	h_1->Divide(h_2);
	h_1->Rebin(10);
	double data_A_normalization_new = h_1->Integral("width");
	h_1->Scale(1/data_A_normalization_new);
	TCanvas *c1 = new TCanvas("","",800,600);
	c1->Divide(2,2);
	c1->cd(1);
	h_1->GetXaxis()->SetRangeUser(0,0.4);
	h_1->Draw("HIST");
	c1->cd(2);
	h_3->GetXaxis()->SetRangeUser(0,0.4);
	h_3->Draw();
	c1->cd(3);
	h_2->GetXaxis()->SetRangeUser(0,0.4);
	h_2->Draw();

	c1->SaveAs("./A_ratio.pdf");
}