void checkhis(int opt = 3){
	//opt: 1 == 0 - 100, 2 = 0 - 24, 3 = 24 - 100
	TFile *file = new TFile("./mc_bk_plot/W_histogram.root","UPDATE");
	file->cd();
	TH1D* h_1;
	TH1D* h_2;
	TH1D* h_3;
	TH1D* h_4;
	TH1D* h_5;
	TH1D* h_6;
	TH1D* h_7;

	TH1D* mc_A;
	TH1D* data_A;

	TH1D* h_DY_weight;
	TH1D* h_tt_weight;
	TH1D* h_W_weight;

	h_DY_weight = (TH1D*) file->Get("DY_weight_factor");
	h_tt_weight = (TH1D*) file->Get("tt_weight_factor");
	h_W_weight = (TH1D*) file->Get("W_weight_factor");

	double w_norm = h_W_weight->Integral();
	double tt_norm = h_tt_weight->Integral();
	double dy_norm = h_DY_weight->Integral();

	double ttbar_XS = 69.0;
    double Wjet_XS = 21159;
    double DY_XS = 2010;

	double muLumi = 1682.8;
    double eLumi = 1670.7 ;
    double netLumi = 1695.6;

	double Nmb = 11775759052;

	float crossSectionModifier = 0.92623216;

	if (opt == 1){
		h_1 = (TH1D*)file->Get("mc_W_bk");
		h_2 = (TH1D*)file->Get("ZTOUU");
		h_3 = (TH1D*)file->Get("mc_bk_tt");
		h_4 = (TH1D*)file->Get("Tau_BK");
		h_5 = (TH1D*)file->Get("data_hist");
		h_6 = (TH1D*)file->Get("data_same_hist");
		h_7 = (TH1D*)file->Get("mc_estimate_ss");
		mc_A = (TH1D*)file->Get("MC_A_ALL");
		data_A = (TH1D*)file->Get("data_A_ALL");
	}
	if (opt == 2){
		h_1 = (TH1D*)file->Get("mc_W_bk_24");
		h_2 = (TH1D*)file->Get("ZTOUU_24");
		h_3 = (TH1D*)file->Get("mc_bk_tt_24");
		h_4 = (TH1D*)file->Get("Tau_BK_24");
		h_5 = (TH1D*)file->Get("data_hist_24");
		h_6 = (TH1D*)file->Get("data_same_hist_24");
		h_7 = (TH1D*)file->Get("mc_estimate_ss_24");
		mc_A = (TH1D*)file->Get("MC_A_24");
		data_A = (TH1D*)file->Get("data_A_24");
	}
	if (opt == 3){
		h_1 = (TH1D*)file->Get("mc_W_bk_100");
		h_2 = (TH1D*)file->Get("ZTOUU_100");
		h_3 = (TH1D*)file->Get("mc_bk_tt_100");
		h_4 = (TH1D*)file->Get("Tau_BK_100");
		h_5 = (TH1D*)file->Get("data_hist_100");
		h_6 = (TH1D*)file->Get("data_same_hist_100");
		h_7 = (TH1D*)file->Get("mc_estimate_ss_100");
		mc_A = (TH1D*)file->Get("MC_A_100");
		data_A = (TH1D*)file->Get("data_A_100");
	}
	if (h_1->GetNbinsX()!= 120 ||h_2->GetNbinsX()!= 120 ||h_3->GetNbinsX()!= 120 ||h_4->GetNbinsX()!= 120 ||h_5->GetNbinsX()!= 120 ||h_6->GetNbinsX()!= 120 ||h_7->GetNbinsX()!= 120 ) cout << "Binning Problems!!! not 120 bins" << endl;

	cout << "The integral of data A is " << data_A->Integral() << endl;
	cout << "The integral of data hist is " << h_5->Integral() << endl;
	cout << "The underflow content is " << data_A->GetBinContent(0) << endl;
	cout << "The overflow content is " << data_A->GetBinContent(1001) << endl;
	cout << "The integral of MC A is " << mc_A->Integral() << endl;
	cout << "The underflow content is " << mc_A->GetBinContent(0) << endl;
	cout << "The overflow content is " << mc_A->GetBinContent(1001) << endl;

    double mc_A_normalization_factor = mc_A->Integral("width");
	double data_A_normalization_factor = data_A->Integral("width");

	mc_A->Scale(1/mc_A_normalization_factor);
	data_A->Scale(1/data_A_normalization_factor);

	h_1->Scale(Wjet_XS*crossSectionModifier/w_norm);
	h_2->Scale(DY_XS*crossSectionModifier/dy_norm);
	h_3->Scale(ttbar_XS*crossSectionModifier/tt_norm);
	h_4->Scale(DY_XS*crossSectionModifier/dy_norm);
	h_7->Scale(DY_XS*crossSectionModifier/dy_norm);

	TH1D* h_data_normalized;
	h_data_normalized = (TH1D*) h_5 ->Clone();

	TH1D* h_mc_normalized;
	h_mc_normalized = (TH1D*)h_2->Clone();
	h_mc_normalized->Add(h_1,1);
	h_mc_normalized->Add(h_3,1);
	h_mc_normalized->Add(h_4,1);
	h_mc_normalized->Add(h_7,1);

	TH1D* h_mc_bk_normalized;
	h_mc_bk_normalized = (TH1D*)h_3->Clone();
	h_mc_bk_normalized->Add(h_4,1);
	h_mc_bk_normalized->Add(h_7,1);
	h_mc_bk_normalized->Add(h_1,1);

	double normalization_factor_mc_bk = h_mc_bk_normalized->Integral("width");
	h_mc_bk_normalized->Scale(1/normalization_factor_mc_bk);
	double normalization_factor_mc_120bins = h_mc_normalized->Integral("width");
	h_mc_normalized->Scale(1/normalization_factor_mc_120bins);
	double normalization_factor_data = h_data_normalized->Integral(1,120,"width");
	h_data_normalized->Scale(1/normalization_factor_data);


	h_1->Rebin(4);
	h_2->Rebin(4);
	h_3->Rebin(4);
	h_4->Rebin(4);
	h_5->Rebin(4);
	h_6->Rebin(4);
	h_7->Rebin(4);
	//h_5->Scale(1/1682.8);
	//h_6->Scale(1/1682.8);

	TH1D* total_mc = new TH1D("total_mc","total_mc",30,60,120);
	total_mc->Add(h_1,1);
	total_mc->Add(h_2,1);
	total_mc->Add(h_3,1);
	total_mc->Add(h_4,1);
	total_mc->Add(h_7,1);

	double normalization_factor_mc = total_mc->Integral("width");
	h_1->Scale(1/normalization_factor_mc);
	h_2->Scale(1/normalization_factor_mc);
	h_3->Scale(1/normalization_factor_mc);
	h_4->Scale(1/normalization_factor_mc);
	h_7->Scale(1/normalization_factor_mc);
	double normalization_factor_data_30bins = h_5->Integral("width");
	h_5->Scale(1/normalization_factor_data_30bins);

	TH1D* total_mc_normalized_30bins = new TH1D("","",30,60,120);
	total_mc_normalized_30bins->Add(h_1);
	total_mc_normalized_30bins->Add(h_2);
	total_mc_normalized_30bins->Add(h_3);
	total_mc_normalized_30bins->Add(h_4);
	total_mc_normalized_30bins->Add(h_7);

	cout << "The area for normalized bk, totalmc, data, rebinned total mc, rebinned data is " << h_mc_bk_normalized->Integral("width") << " " << h_mc_normalized->Integral("width") << " " << h_data_normalized->Integral("width") << " " << total_mc_normalized_30bins->Integral("width") << " " << h_5->Integral("width") << endl; 
	/*h_1->Scale(muLumi/netLumi * Nmb / 0.9 * 5.649e-9* 67.6/70.0);
	h_2->Scale(muLumi/netLumi * Nmb / 0.9 * 5.649e-9 * 67.6/70.0);
	h_3->Scale(muLumi/netLumi * Nmb / 0.9 * 5.649e-9 * 67.6/70.0);
	h_4->Scale(muLumi/netLumi * Nmb / 0.9 * 5.649e-9 * 67.6/70.0);
	h_7->Scale(muLumi/netLumi * Nmb / 0.9 * 5.649e-9 * 67.6/70.0);*/

	cout << "W is " << Wjet_XS*crossSectionModifier/w_norm << " DY is " << DY_XS*crossSectionModifier/dy_norm << " TT is " << ttbar_XS*crossSectionModifier/tt_norm << endl; 
	//h_5->Scale();
	//h_6->Scale();

	TCanvas *c1 = new TCanvas("","",800,600);
	c1->SetLogy();
	c1->cd();
	h_2->SetFillColor(kOrange+1);
	h_1->SetFillColor(kMagenta+2);
	h_3->SetFillColor(kGray);
	h_4->SetFillColor(kGreen);
	h_7->SetFillColor(kBlue);
	h_1->SetOption("HIST");
	h_2->SetOption("HIST");
	h_3->SetOption("HIST");
	h_4->SetOption("HIST");
	h_7->SetOption("HIST");
	h_5->SetMarkerStyle(20); // Marker style (circle)
    h_5->SetMarkerSize(0.5); // Marker size
    h_5->SetMarkerColor(kRed); 

	TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Specify the legend position
	legend->AddEntry(h_5, "data_18", "p"); // Add an entry for histo2
	legend->AddEntry(h_2, "MC_Z->u+u-", "f"); // Add an entry for histo2
	legend->AddEntry(h_7,"Same_sign(QCD)","f");
	legend->AddEntry(h_4, "MC_Z->tau+tau-", "f"); // Add an entry for histo2
    legend->AddEntry(h_1, "MC_W_background", "f"); // Add an entry for histo1
	legend->AddEntry(h_3, "MC_ttbar_background", "f"); // Add an entry for histo2
    //legend->SetBorderSize(0); // Set border size of the legend

	THStack *s1 = new THStack("s1","");
	s1->SetTitle("Z_Mass_Spectrum");
	s1->SetMinimum(1e-5);
	s1->SetMaximum(2e-1);
	//s1->Add(h_7,"HIST");
	s1->Add(h_3,"HIST");
	s1->Add(h_1,"HIST");
	s1->Add(h_4,"HIST");
	s1->Add(h_7,"HIST");
	s1->Add(h_2,"HIST");
	s1->Draw();
	s1->GetXaxis()->SetTitle("Invariant mass (GeV)");
	s1->GetYaxis()->SetTitle("Counts/4.0");
	h_5->Draw("SAME");
	legend->Draw("SAME");
	if (opt == 1)c1->SaveAs("./compositeplot/Composite.pdf");
	if (opt == 2)c1->SaveAs("./compositeplot/Composite_24.pdf");
	if (opt == 3)c1->SaveAs("./compositeplot/Composite_100.pdf");


	TCanvas *c2 = new TCanvas("","",800,600);
	c2->cd();
	mc_A->SetLineColor(kBlue);
	data_A->SetLineColor(kRed);
	mc_A->GetXaxis()->SetRangeUser(0,0.1);
	data_A->GetXaxis()->SetRangeUser(0,0.1);
	data_A->Draw("P");
	mc_A->Draw("P SAME");

	if (opt == 1) c2->SaveAs("./compositeplot/A_all.pdf");
	if (opt == 2) c2->SaveAs("./compositeplot/A_24.pdf");
	if (opt == 3) c2->SaveAs("./compositeplot/A_100.pdf");

	file->cd();
	if (opt == 1){
		h_mc_normalized->Write("Normalized_mc",2);
		h_data_normalized->Write("Normalized_data",2);
		h_mc_bk_normalized->Write("Normalized_mc_bk",2);
	}
	if (opt == 2){
		h_mc_normalized->Write("Normalized_mc_24",2);
		h_data_normalized->Write("Normalized_data_24",2);
		h_mc_bk_normalized->Write("Normalized_mc_bk_24",2);
	}
	if (opt == 3){
		h_mc_normalized->Write("Normalized_mc_100",2);
		h_data_normalized->Write("Normalized_data_100",2);
		h_mc_bk_normalized->Write("Normalized_mc_bk_100",2);
	}
	file->Close();

}