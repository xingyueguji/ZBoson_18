class plotting_helper{
	public:

	plotting_helper();
	~plotting_helper();
	void areanormalize(TH1D* h_1);
	void luminormalize(TH1D* h_1,int opt, double weight);
	void compositeplot(TH1D* h_1, TH1D* h_2, TH1D* h_3, TH1D* h_4, TH1D* h_5, TH1D* h_6,int x);
	void acoplot(TH1D* h_1, TH1D* h_2, int x);
	void savehistogram(TH1D* h_1, TH1D* h_2, TH1D* h_3, int x, TFile *f1);

	double ttbar_XS = 69.0;
    double Wjet_XS = 21159;
    double DY_XS = 2010;

	double muLumi = 1682.8;
    double eLumi = 1670.7 ;
    double netLumi = 1695.6;

	double Nmb = 11775759052;

	float crossSectionModifier = 0.92623216;


};


plotting_helper::plotting_helper(){
}

void plotting_helper::areanormalize(TH1D* h_1){
	double normalization_factor = h_1->Integral("width");
	h_1->Scale(1/normalization_factor);
}

void plotting_helper::luminormalize(TH1D* h_1,int opt,double weight){
	//opt 1 == DY, opt 2 == W, opt 3 = TT
	if (opt == 1) h_1->Scale(DY_XS*crossSectionModifier/weight);
	if (opt == 2) h_1->Scale(Wjet_XS*crossSectionModifier/weight);
	if (opt == 3) h_1->Scale(ttbar_XS*crossSectionModifier/weight);
}

void plotting_helper::compositeplot(TH1D* h_1, TH1D* h_2, TH1D* h_3, TH1D* h_4, TH1D* h_5, TH1D* h_6,int x){
	TCanvas *c1 = new TCanvas("","",800,600);
	c1->SetLogy();
	c1->cd();
	h_1->SetMarkerStyle(20); // Marker style (circle)
    h_1->SetMarkerSize(0.5); // Marker size

    h_1->SetMarkerColor(kRed); 
	h_2->SetFillColor(kOrange+1);
	h_3->SetFillColor(kBlue);
	h_4->SetFillColor(kGreen);
	h_5->SetFillColor(kMagenta+2);
	h_6->SetFillColor(kGray);

	h_2->SetOption("HIST");
	h_3->SetOption("HIST");
	h_4->SetOption("HIST");
	h_5->SetOption("HIST");
	h_6->SetOption("HIST");

	TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); 
	legend->AddEntry(h_1, "data_18 - EM background", "p"); 
	legend->AddEntry(h_2, "MC_Z->u+u-", "f"); 
	legend->AddEntry(h_3,"Same_sign(QCD)","f");
	legend->AddEntry(h_4, "MC_Z->tau+tau-", "f"); 
    legend->AddEntry(h_5, "MC_W_background", "f");
	legend->AddEntry(h_6, "MC_ttbar_background", "f"); 
    //legend->SetBorderSize(0); // Set border size of the legend
	THStack *s1 = new THStack("s1","");
	s1->SetTitle("Z_Mass_Spectrum");
	s1->SetMinimum(1e-5);
	s1->SetMaximum(2e-1);
	//s1->Add(h_7,"HIST");
	s1->Add(h_6,"HIST");
	s1->Add(h_5,"HIST");
	s1->Add(h_4,"HIST");
	s1->Add(h_3,"HIST");
	s1->Add(h_2,"HIST");
	s1->Draw();
	s1->GetXaxis()->SetTitle("Invariant mass (GeV)");
	s1->GetYaxis()->SetTitle("Counts/4.0");
	h_1->Draw("SAME");
	legend->Draw("SAME");
	c1->SaveAs(Form("./newcomposite/composite_%i.pdf",x));

}
void plotting_helper::acoplot(TH1D* h_1, TH1D* h_2, int x){
	TCanvas *c2 = new TCanvas("","",800,600);
	c2->cd();
	h_1->SetLineColor(kBlue);
	h_2->SetLineColor(kRed);
	h_1->GetXaxis()->SetRangeUser(0,0.1);
	h_2->GetXaxis()->SetRangeUser(0,0.1);
	h_2->Draw("P");      // data
	h_1->Draw("P SAME"); // mc

    c2->SaveAs(Form("./acoplot/A_%i.pdf",x));
}

void plotting_helper::savehistogram(TH1D* h_1, TH1D* h_2, TH1D* h_3, int x, TFile *f1){
	f1->cd();
	h_1->Write(Form("Normalized_mc_%i",x),2);
	h_2->Write(Form("Normalized_data_%i",x),2);
	h_3->Write(Form("Normalized_mc_bk_%i",x),2);

}



