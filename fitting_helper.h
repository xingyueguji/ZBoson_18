class fitting_helper{
	public:

	fitting_helper();
	~fitting_helper();
	void gethistograms(int opt);
	void setupfittingfunction(int cent);
	void fitting(int opt1);
	void plotting(int opt1, int opt, int cent);
	void setupfittingfunctionchi();
	void plotTGraph(TGraphErrors* x1, TGraphErrors* x2, TGraphErrors *x3, int opt, int opt1);

	TString raw_file = "./rootfile/normalized/rawfile.root";
	TString ycut_file = "./rootfile/normalized/ycutfile.root";
	TString eff_file = "./rootfile/normalized/efffile.root";
	TString ycut_eff_file = "./rootfile/normalized/ycut_eff_file.root";
	TString etacut_file = "./rootfile/normalized/etacut_file.root";
	TString etacut_eff_file = "./rootfile/normalized/etacut_eff_file.root";

	TH1D* h_normalized_mc[11];
	TH1D* h_normalized_data[11];
	TH1D* h_normalized_mc_bk[11];

	TFile* f1;
	TFile* f2;

	double bwmeanip = 90.8513, widthip = 2.911, cbmeanip = 0, sigmaip = 1.1104, alphaip = 1.828, nip = 0.972, decayip = -0.1163, fsigip = 0.99199, sig1fracip = 0.99199;

	//This is Breit_Wigner
	RooRealVar *x;
	RooRealVar *conwindow;
	RooRealVar *bwmean;
	//RooRealVar m0("m0", "m0",m0ip,80,100);
	RooRealVar *width;
	RooBreitWigner *bw;
	//This is CrystalBall
	RooRealVar *cbmean;
	RooRealVar *cbsigma;
	RooRealVar *cbalpha;
	RooRealVar *cbn;
	RooCBShape *cb;
	//This is Exponential
	RooRealVar *decay;
	RooExponential *exp;
	//This is BW+CB
	RooNumConvPdf *pdf;
	//This is BW+CB+EXP
	RooRealVar *fsig;
	RooAddPdf *purepdf;
	//This is BW+CB+Template
	RooDataHist *background;
	RooHistPdf *histpdf1;
	RooRealVar *sig1frac;
	RooAddPdf *temp_sig;
	RooDataHist *data;
	RooDataHist *mc;
	TCanvas *c1;
	RooPlot* pullFrame;
	RooHist* residuals;
	TPaveText *textBox;
	RooPlot* frame;
	RooPlot* framecheck;

	RooFFTConvPdf* newconvpdf;

	RooRealVar *convocenter;
	RooRealVar *convowidth;

	double meanshift = 0;
	double meanshifterror = 0;

	double widthshift = 0;
	double widtherror = 0;

	Int_t h_normalized_data_entry[11] = {};







};

fitting_helper::fitting_helper(){
}

void fitting_helper::gethistograms(int opt){
	if (opt == 1)f1 = new TFile(raw_file, "READ");
	if (opt == 2)f1 = new TFile(ycut_file, "READ");
	if (opt == 3)f1 = new TFile(eff_file, "READ");
	if (opt == 4)f1 = new TFile(ycut_eff_file, "READ");
	if (opt == 5)f1 = new TFile(etacut_file,"READ");
	if (opt == 6)f1 = new TFile(etacut_eff_file,"READ");
	//f2 = new TFile(mc_file_path,"READ");

	for (int i=0; i<11; i++){
		h_normalized_mc[i] = (TH1D*) f1->Get(Form("Normalized_mc_%i",i));
		h_normalized_data[i] = (TH1D*) f1->Get(Form("Normalized_data_%i",i));
		h_normalized_mc_bk[i] = (TH1D*) f1->Get(Form("Normalized_mc_bk_%i",i));
		h_normalized_data_entry[i] = h_normalized_data[i]->GetEntries();
	}

	
}

void fitting_helper::setupfittingfunction(int cent){
	convocenter = new RooRealVar("convocenter", "", 91);
	convowidth = new RooRealVar("convowidth","convowidth",500);

	x = new RooRealVar("x", "x",60., 120. );
	x->setBinning(RooBinning(10000,60,120),"cache");
	conwindow = new RooRealVar("conwindow","conwindow",-10,125);
	bwmean = new RooRealVar("bwmean", "bwmean",bwmeanip,80,100);
	//RooRealVar m0("m0", "m0",m0ip,80,100);
	width = new RooRealVar("width", "width",widthip,0,10 );
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width );
	cbmean = new RooRealVar("mean","mean",cbmeanip);
	cbsigma = new RooRealVar("sigma","sigma",sigmaip,0,10);
	cbalpha = new RooRealVar("alpha","alpha",alphaip,-5,5);
	cbn = new RooRealVar("n","n", nip ,-5,5);
	cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	decay = new RooRealVar("decay","decay",decayip,-5,+5);
	exp = new RooExponential("exp", "Exponential PDF", *x, *decay);
	newconvpdf = new RooFFTConvPdf("newconvpdf","newconvpdf", *x, *bw, *cb);
	//pdf = new RooNumConvPdf("pdf","pdf", *x, *bw, *cb);
	fsig = new RooRealVar("fsig","signal fraction",fsigip,0,1);
    purepdf = new RooAddPdf("purepdf","cbconbw+exp",RooArgList(*newconvpdf,*exp),*fsig);
	background = new RooDataHist("background","background",*x,h_normalized_mc_bk[cent]);
	histpdf1 = new RooHistPdf("histpdf1", "histpdf1", *x, *background, 0);
	sig1frac = new RooRealVar("sig1frac", "fraction of component 1 in signal", sig1fracip,0,1);
	temp_sig = new RooAddPdf("temp_sig","temp_sig",RooArgList(*newconvpdf,*histpdf1),*sig1frac);
	data = new RooDataHist("data","data",*x,h_normalized_data[cent]);
	mc = new RooDataHist("mc","mc",*x,h_normalized_mc[cent]);

	//pdf->setConvolutionWindow(*convocenter, *convowidth, 1);
	//RooFit::Minimizer("Minuit","migrad")
}

void fitting_helper::fitting(int opt1){
	//pdf->fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	//if (opt1 == 1 )pdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true));
	if (opt1 == 1 )newconvpdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
	//purepdf->fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	if (opt1 == 2 )purepdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
	//temp_sig->fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	if (opt1 == 3)temp_sig->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
}
void fitting_helper::plotting(int opt1, int opt, int cent){
	frame = x->frame();
	framecheck = x->frame();
	if (opt == 1) {
		data->plotOn(frame,RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		newconvpdf->plotOn(frame,RooFit::LineWidth(1));
		newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
	}
	if (opt == 2) {
		data->plotOn(frame,RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		purepdf->plotOn(frame,RooFit::LineWidth(1));
		purepdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		purepdf->plotOn(framecheck,RooFit::Components(*exp), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	}
	if (opt == 3) {
		data->plotOn(frame,RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		temp_sig->plotOn(frame,RooFit::LineWidth(1));
		temp_sig->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		temp_sig->plotOn(framecheck,RooFit::Components(*histpdf1), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	}
	residuals = frame->residHist("","",true,false); //true = pull, false = center of the bin
	pullFrame = x->frame();
	pullFrame->addPlotable(residuals, "P");
	pullFrame->SetTitle("");
	pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
	frame->SetTitle("");
	frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
	c1 = new TCanvas("c1","",1600,1200);
	c1->Divide(2,2);
	c1->cd(1);
	frame->Draw();
	double meanofBW = bwmean->getVal() - 91.1876;
	double meanerrorofBW = bwmean->getError();
	double widthofBW = width->getVal() - 2.4955;
	double widtherrorofBW = width->getError();
	Int_t numEntries = data->numEntries();
	
	meanshift = meanofBW;
	meanshifterror = meanerrorofBW;
	widthshift = widthofBW;
	widtherror = widtherrorofBW;

	textBox = new TPaveText(0.1, 0.5, 0.3, 0.8, "NDC");
    textBox->SetFillColor(0);
    textBox->SetTextSize(0.03);
    textBox->AddText(Form("Mean value: %.2f", meanofBW));
	textBox->AddText(Form("Width value: %.2f",widthofBW));
	textBox->AddText(Form("# of Entries: %i",h_normalized_data_entry[cent]));
    textBox->Draw();
	c1->cd(2);
	pullFrame->Draw();
	c1->cd(3);
	framecheck->Draw();
	if (opt1 == 1){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/raw/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/raw/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/raw/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 2){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/ycut/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/ycut/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/ycut/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 3){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 4){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/ycut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/ycut_eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/ycut_eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 5){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/etacut/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/etacut/cb+bw_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/etacut/cb+bw_%i.pdf",cent));

	}
	if (opt1 == 6){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/etacut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/etacut_eff/cb+bw_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/etacut_eff/cb+bw_%i.pdf",cent));
	}




}

void fitting_helper::setupfittingfunctionchi(){
	meanshift = 2;
}

void fitting_helper::plotTGraph(TGraphErrors* x1, TGraphErrors* x2, TGraphErrors *x3, int opt, int opt1){
	gStyle->SetEndErrorSize(6);
	gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.11);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
	gStyle->SetTitleFont(42);
	gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.03, "XYZ");

	TCanvas* c_2 = new TCanvas("","",600,600);
	TPad *pad = new TPad("pad", "Pad", 0, 0, 1, 1);
    pad->Draw();
    pad->cd();
	x1->SetTitle("");
	x1->GetXaxis()->SetTitle("Central centrality");
	x1->GetXaxis()->CenterTitle();
	if (opt == 1) 	x1->GetYaxis()->SetTitle("dM(GeV)");
	if (opt == 2) 	x1->GetYaxis()->SetTitle("dWidth(GeV)");
	x1->GetYaxis()->CenterTitle();
	x1->GetYaxis()->SetTitleOffset(1.9);
	
	if (opt == 1) x1->GetYaxis()->SetRangeUser(-0.6,0);
	if (opt == 2) x1->GetYaxis()->SetRangeUser(-0.2,1.4);

	x1->SetLineWidth(2);
	x2->SetLineWidth(2);
	x3->SetLineWidth(2);

	x1->SetMarkerSize(1.7);
	x1->SetMarkerStyle(25);
	x2->SetMarkerSize(1.7);
	x2->SetMarkerStyle(4);
	x3->SetMarkerSize(1.7);
	x3->SetMarkerStyle(kOpenTriangleUp);

	x1->SetMarkerColor(2);
	x1->SetLineColor(2);
	x2->SetMarkerColor(4);
	x2->SetLineColor(4);
	x3->SetMarkerColor(8);
	x3->SetLineColor(8);


	TF1 *f1 = new TF1("f1", "[0]", 0, 80);
	TF1 *f2 = new TF1("f2", "[0]", 0, 80);
	TF1 *f3 = new TF1("f3", "[0]", 0, 80);

	f1->SetLineColor(2);
	f2->SetLineColor(4);
	f3->SetLineColor(8);

	TLegend *legend = new TLegend(0.6, 0.6, 0.8, 0.9);
	legend->SetTextFont(42);
	legend->SetBorderSize(0);
    legend->AddEntry(x1, "cb+bw", "lep"); 
	legend->AddEntry(x2, "cb+bw+exp", "lep"); 
	legend->AddEntry(x3, "cb+bw+template", "lep"); 
	legend->AddEntry(f1,"cb+bw fitted","l");
	legend->AddEntry(f2,"cb+bw+exp fitted","l");
	legend->AddEntry(f3,"cb+bw+template fitted","l");
    legend->SetFillColor(kWhite);



	x1->Draw("AP");
	x2->Draw("P");
	x3->Draw("P");
	x1->Fit(f1,"QR");
	x2->Fit(f2,"QR");
	x3->Fit(f3,"QR");

	double fittedValue = f1->GetParameter(0);
	double fittedValue2 = f2->GetParameter(0);
	double fittedValue3 = f3->GetParameter(0);

	TPaveText *pt = new TPaveText(0.2, 0.8, 0.55, 0.9, "brNDC"); // normalized coordinates
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextAlign(12); // Align left and vertically centered
	pt->SetTextFont(42);
	pt->SetTextSize(0.02);
	pt->SetMargin(0.02);
	pt->AddText(Form("Fitted constant for cb+bw = %.4f", fittedValue));
	pt->AddText(Form("Fitted constant for cb+bw+exp = %.4f", fittedValue2));
	pt->AddText(Form("Fitted constant for cb+bw+temp = %.4f", fittedValue3));

	legend->Draw();
	pt->Draw();

	if (opt == 1 && opt1 == 1)	c_2->SaveAs("./dM_dSig/raw/dM.pdf");
	if (opt == 2 && opt1 == 1) 	c_2->SaveAs("./dM_dSig/raw/dSig.pdf");

	if (opt == 1 && opt1 == 2)	c_2->SaveAs("./dM_dSig/ycut/dM.pdf");
	if (opt == 2 && opt1 == 2) 	c_2->SaveAs("./dM_dSig/ycut/dSig.pdf");

	if (opt == 1 && opt1 == 3)	c_2->SaveAs("./dM_dSig/eff/dM.pdf");
	if (opt == 2 && opt1 == 3) 	c_2->SaveAs("./dM_dSig/eff/dSig.pdf");

	if (opt == 1 && opt1 == 4)	c_2->SaveAs("./dM_dSig/ycut_eff/dM.pdf");
	if (opt == 2 && opt1 == 4) 	c_2->SaveAs("./dM_dSig/ycut_eff/dSig.pdf");

	if (opt == 1 && opt1 == 5)	c_2->SaveAs("./dM_dSig/etacut/dM.pdf");
	if (opt == 2 && opt1 == 5) 	c_2->SaveAs("./dM_dSig/etacut/dSig.pdf");

	if (opt == 1 && opt1 == 6)	c_2->SaveAs("./dM_dSig/etacut_eff/dM.pdf");
	if (opt == 2 && opt1 == 6) 	c_2->SaveAs("./dM_dSig/etacut_eff/dSig.pdf");

}
