void checkmcfitsingle(int cent = 2){
	TFile *f1;
	TString mcfile = "./rootfile/mc_signal.root";

	f1 = new TFile(mcfile,"READ");

	RooDataSet* rooreco[11];
	RooDataSet* rooreco_y[11];
	RooDataSet* rooreco_eff[11];
	RooDataSet* rooreco_y_eff[11];
	RooDataSet* rooreco_eta[11];
	RooDataSet* rooreco_eta_eff[11];

	RooDataSet* roogen[11];
	RooDataSet* roogen_y[11];
	RooDataSet* roogen_eta[11];

	//This is for data //double bwmeanip = 90.8513, widthip = 2.911, cbmeanip = 0, sigmaip = 1.1104, alphaip = 1.828, nip = 0.972, decayip = -0.1163, fsigip = 0.99199, sig1fracip = 0.99199;
	double bwmeanip = 91.0313, widthip = 2.6933, cbmeanip = 0, sigmaip = 1.1152, alphaip = 1.7402, nip = 1.0894, decayip = -0.1163, fsigip = 0.99199, sig1fracip = 0.99199;
	// This is MC above.
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
	RooFitResult* fitresult;

	RooFFTConvPdf* newconvpdf;

	RooRealVar *convocenter;
	RooRealVar *convowidth;

	RooDataSet *unbinned;
	RooArgSet *useforunbinned;

	RooRealVar *weight;

	Int_t numEntries;
	Double_t integral;

	for (int i=0; i<11; i++){
		rooreco[i] = (RooDataSet*)f1->Get(Form("rooreco_%i",i));
		rooreco_y[i] = (RooDataSet*)f1->Get(Form("rooreco_y_%i",i));
		rooreco_eff[i] = (RooDataSet*)f1->Get(Form("rooreco_eff_%i",i));
		rooreco_y_eff[i] = (RooDataSet*)f1->Get(Form("rooreco_y_eff_%i",i));
		rooreco_eta[i] = (RooDataSet*)f1->Get(Form("rooreco_eta_%i",i));
		rooreco_eta_eff[i] = (RooDataSet*)f1->Get(Form("rooreco_eta_eff_%i",i));

		roogen[i] = (RooDataSet*)f1->Get(Form("roogen_%i",i));
		roogen_y[i] = (RooDataSet*)f1->Get(Form("roogen_y_%i",i));
		roogen_eta[i] = (RooDataSet*)f1->Get(Form("roogen_eta_%i",i));
	}

	x = new RooRealVar("roomass","roomass",60,120);
	x->setBinning(RooBinning(10000,60,120),"cache");
	bwmean = new RooRealVar("bwmean", "bwmean",bwmeanip,80,100);
	width = new RooRealVar("width", "width",widthip,0,10 );
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width );
	cbmean = new RooRealVar("mean","mean",cbmeanip);
	cbsigma = new RooRealVar("sigma","sigma",sigmaip,0,10);
	cbalpha = new RooRealVar("alpha","alpha",alphaip,-5,5);
	cbn = new RooRealVar("n","n", nip ,-5,5);
	cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	//decay = new RooRealVar("decay","decay",decayip,-5,+5);
	//exp = new RooExponential("exp", "Exponential PDF", *x, *decay);
	newconvpdf = new RooFFTConvPdf("newconvpdf","newconvpdf", *x, *bw, *cb); // This is CB + BW
	//fsig = new RooRealVar("fsig","signal fraction",fsigip,0,1);
    //purepdf = new RooAddPdf("purepdf","cbconbw+exp",RooArgList(*newconvpdf,*exp),*fsig); // This is EXP
	//background = new RooDataHist("background","background",*x,h_normalized_mc_bk[cent]);
	//histpdf1 = new RooHistPdf("histpdf1", "histpdf1", *x, *background, 0);
	//sig1frac = new RooRealVar("sig1frac", "fraction of component 1 in signal", sig1fracip,0,1);
	//temp_sig = new RooAddPdf("temp_sig","temp_sig",RooArgList(*newconvpdf,*histpdf1),*sig1frac); // This is template

	fitresult = newconvpdf->fitTo(*rooreco_eff[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
	frame = x->frame(RooFit::Title("Z mass fit"));

	rooreco_eff[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
	newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
	newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
	numEntries = rooreco_eff[cent]->numEntries();
	integral = rooreco_eff[cent]->sumEntries();
	residuals = frame->residHist("roodata","fit",true,false);
	pullFrame = x->frame();
	pullFrame->addPlotable(residuals, "P");
	pullFrame->SetTitle("");
	pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
	frame->SetTitle("");
	frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
	c1 = new TCanvas("c1","",2200,800);
	c1->Divide(2,1);
	c1->cd(1);
	frame->Draw();
	c1->cd(2);
	pullFrame->Draw();

	c1->SaveAs(Form("./singlemcfitcheck/rooreco_eff_%i.pdf",cent));





}