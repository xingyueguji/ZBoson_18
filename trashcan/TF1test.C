void TF1test(){
	/*TF1* f_1 = new TF1("f_1","[0]*ROOT::Math::crystalball_function(x - [4], [1], [2], [3])",60,120); // constant, alpha,N,sig,mean
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
	c1->SaveAs("f1_test.pdf");*/

RooRealVar t("t","t",-10,30) ;
// Theoretical model
RooRealVar ml("ml","mean landau",5.,-20,20) ;
RooRealVar sl("sl","sigma landau",1,0.1,10) ;
RooLandau landau("lx","lx",t,ml,sl) ;
// Detector response function
RooRealVar mg("mg","mg",0) ;
RooRealVar sg("sg","sg",2,0.1,10) ;
RooGaussian gauss("gauss","gauss",t,mg,sg) ;
// Define sampling frequency
//t.setBins("fft",10000) ;
// Construct convolution
RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;
// Sample 1000 events in x from gxlx
RooDataSet* data = lxg.generate(t,10000) ;
// Fit gxlx to data
lxg.fitTo(*data) ;
// Plot data, fitted p.d.f
//RooPlot* frame = t.frame(RooFit::Title("landau (x) gauss convolution")) ;
//data->plotOn(frame) ;
//lxg.plotOn(frame) ;
//landau.plotOn(frame,RooFit::LineStyle(kDashed)) ;

}