void TF1test(){
	TF1* f_1 = new TF1("f_1","[0]*ROOT::Math::crystalball_function(x - [4], [1], [2], [3])",60,120); // constant, alpha,N,sig,mean
	TF1* f_2 = new TF1("f_2","[0]*ROOT::Math::breitwigner_pdf (x-91.1876-[1],[2]+2.4955)",60,120); // constant, mean, width
	f_1->SetParameters(1,1,1,5,91);

	TF1Convolution *f_conv = new TF1Convolution("f_1", "f_2", -5, 120);
	TF1 *f = new TF1("f", *f_conv, 60., 120., f_conv->GetNpar());
	f_conv->SetRange(60., 120.);
	//f->SetParameters(1., -0.3, 0., 1.);
	f->FixParameter(4,0);
	f->SetParameters(1,1.8,0.9378,1,0,0.69157,-0.474819,1.14813);

	TFile *f1 = new TFile("./rootfile/data.root","READ");
	TH1D* h_1;
	h_1 = (TH1D*) f1->Get("Normalized_data_1");
	h_1->Fit(f);

	TCanvas *c1 = new TCanvas("","",800,600);
	c1->cd();
	h_1->Draw();
	c1->SaveAs("f1_test.pdf");

}