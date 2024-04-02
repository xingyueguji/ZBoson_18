class fitting_helper{
	public:

	fitting_helper();
	~fitting_helper();
	void gethistograms();
	void setupfittingfunction(int cent);
	void fitting(int opt1);
	void plotting(int opt, int cent);
	void setupfittingfunctionchi();

	TString mc_file_path = "./rootfile/mc_bk.root";
	TString data_file_path = "./rootfile/data.root";

	TH1D* h_normalized_mc[11];
	TH1D* h_normalized_data[11];
	TH1D* h_normalized_mc_bk[11];

	TFile* f1;
	TFile* f2;

	double bwmeanip = 90.892, widthip = 3.43133, cbmeanip = 0, sigmaip = 0.544, alphaip = 1.828, nip = 0.9378, decayip = -0.07613, fsigip = 0.99, sig1fracip = 0.99;

	//This is Breit_Wigner
	RooRealVar *x;
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

	double meanshift = 0;
	double meanshifterror = 0;







};

fitting_helper::fitting_helper(){
}

void fitting_helper::gethistograms(){
	f1 = new TFile(data_file_path, "READ");
	//f2 = new TFile(mc_file_path,"READ");

	for (int i=0; i<11; i++){
		h_normalized_mc[i] = (TH1D*) f1->Get(Form("Normalized_mc_%i",i));
		h_normalized_data[i] = (TH1D*) f1->Get(Form("Normalized_data_%i",i));
		h_normalized_mc_bk[i] = (TH1D*) f1->Get(Form("Normalized_mc_bk_%i",i));
	}

	
}

void fitting_helper::setupfittingfunction(int cent){
	x = new RooRealVar("x", "x",60., 120. );
	bwmean = new RooRealVar("bwmean", "bwmean",bwmeanip,80,100);
	//RooRealVar m0("m0", "m0",m0ip,80,100);
	width = new RooRealVar("width", "width",widthip,0,10 );
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width );
	cbmean = new RooRealVar("mean","mean",cbmeanip);
	cbsigma = new RooRealVar("sigma","sigma",sigmaip,0,10);
	cbalpha = new RooRealVar("alpha","alpha",alphaip,-5,5);
	cbn = new RooRealVar("n","n", nip ,-5,5);
	cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	decay = new RooRealVar("decay","decay",decayip,-5,5);
	exp = new RooExponential("exp", "Exponential PDF", *x, *decay);
	pdf = new RooNumConvPdf("pdf","pdf", *x, *bw, *cb);
	fsig = new RooRealVar("fsig","signal fraction",fsigip,0,1);
    purepdf = new RooAddPdf("purepdf","cbconbw+exp",RooArgList(*pdf,*exp),*fsig);
	background = new RooDataHist("background","background",*x,h_normalized_mc_bk[cent]);
	histpdf1 = new RooHistPdf("histpdf1", "histpdf1", *x, *background, 0);
	sig1frac = new RooRealVar("sig1frac", "fraction of component 1 in signal", sig1fracip,0,1);
	temp_sig = new RooAddPdf("temp_sig","temp_sig",RooArgList(*pdf,*histpdf1),*sig1frac);
	data = new RooDataHist("data","data",*x,h_normalized_data[cent]);
	mc = new RooDataHist("mc","mc",*x,h_normalized_mc[cent]);
	//RooFit::Minimizer("Minuit","migrad")
}

void fitting_helper::fitting(int opt1){
	//pdf->fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	if (opt1 == 1 )pdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true));
	//purepdf->fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	if (opt1 == 2 )purepdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true));
	//temp_sig->fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	if (opt1 == 3)temp_sig->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true));
}
void fitting_helper::plotting(int opt, int cent){
	frame = x->frame();
	framecheck = x->frame();
	if (opt == 1) {
		data->plotOn(frame,RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		pdf->plotOn(frame,RooFit::LineWidth(1));
		pdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
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
	c1 = new TCanvas("c1","",1600,1200);
	c1->Divide(2,2);
	c1->cd(1);
	frame->Draw();
	double meanofBW = bwmean->getVal() - 91.1876;
	double meanerrorofBW = bwmean->getError();
	double widthofBW = width->getVal() - 2.4955;
	
	meanshift = meanofBW;
	meanshifterror = meanerrorofBW;

	textBox = new TPaveText(0.1, 0.5, 0.3, 0.8, "NDC");
    textBox->SetFillColor(0);
    textBox->SetTextSize(0.03);
    textBox->AddText(Form("Mean value: %.2f", meanofBW));
	textBox->AddText(Form("Width value: %.2f",widthofBW));
    textBox->Draw();
	c1->cd(2);
	pullFrame->Draw();
	c1->cd(3);
	framecheck->Draw();

	if (opt == 1) c1->SaveAs(Form("./fitresultplot/cb+bw_%i.pdf",cent));
	if (opt == 2) c1->SaveAs(Form("./fitresultplot/cb+bw+exp_%i.pdf",cent));
	if (opt == 3) c1->SaveAs(Form("./fitresultplot/cb+bw+temp_%i.pdf",cent));


}

void fitting_helper::setupfittingfunctionchi(){
	meanshift = 2;
}
