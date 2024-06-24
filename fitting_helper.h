class fitting_helper{
	public:

	fitting_helper();
	~fitting_helper();
	void gethistograms(int opt,bool useunbin);
	void setupfittingfunction(int cent, int opt, bool useunbin);
	RooFitResult* fitting(int opt1,int cent,bool useunbin, bool usew2unbin);
	void plotting(int opt1, int opt, int cent, bool useunbin, bool usew2unbin);
	void setupfittingfunctionchi();
	void plotTGraph(TGraphErrors* x1, TGraphErrors* x2, TGraphErrors *x3, int opt, int opt1,bool useunbin, bool usew2unbin);


	//Below are functions of mc only
	void getmcdataset();
	void setupfittingfunctionmc();
	RooFitResult* fittingmc(int type, bool iseff, bool usew2unbin, int cent);
	void plottingmc(int type, int cent, bool iseff, bool usew2unbin);
	void plotTGraphmc(TGraphErrors* x1, int type, bool iseff, bool usew2unbin, int param, int cent);

	TString raw_file = "./rootfile/normalized/rawfile.root";
	TString ycut_file = "./rootfile/normalized/ycutfile.root";
	TString eff_file = "./rootfile/normalized/efffile.root";
	TString ycut_eff_file = "./rootfile/normalized/ycut_eff_file.root";
	TString etacut_file = "./rootfile/normalized/etacut_file.root";
	TString etacut_eff_file = "./rootfile/normalized/etacut_eff_file.root";
	TString data_file = "./rootfile/data_file.root";
	TString mc_file = "./rootfile/mc_signal.root";

	TH1D* h_normalized_mc[11];
	TH1D* h_normalized_data[11];
	TH1D* h_normalized_mc_bk[11];

    RooDataSet* roodata[11];
	RooDataSet* roodata_y[11];
	RooDataSet* roodata_y_eff[11];
	RooDataSet* roodata_eta[11];
	RooDataSet* roodata_eta_eff[11];
	RooDataSet* roodata_eff[11];

	RooDataSet* rooreco[11];
	RooDataSet* rooreco_y[11];
	RooDataSet* rooreco_eff[11];
	RooDataSet* rooreco_y_eff[11];
	RooDataSet* rooreco_eta[11];
	RooDataSet* rooreco_eta_eff[11];

	RooDataSet* roogen[11];
	RooDataSet* roogen_y[11];
	RooDataSet* roogen_eta[11];

	TFile* f1;
	TFile* f2;
	TFile* f3;

	//This is for data 
	double bwmeanip = 90.8513, widthip = 2.911, cbmeanip = 0, sigmaip = 1.1104, alphaip = 1.828, nip = 0.972, decayip = -0.1163, fsigip = 0.99199, sig1fracip = 0.99199;

	// This is MC.
	double bwmeanrecoip = 91.0313, widthrecoip = 2.6933, cbmeanrecoip = 0, sigmarecoip = 1.2, alpharecoip = 1.7402, nrecoip = 1.0894; 
	double bwmeangenip = 91.32, widthgenip = 2.286, cbmeangenip = 0, sigmagenip = 0.0089,alphagenip = 0.4, ngenip = 1.315; 
	

	RooRealVar *x;
	//This is Breit_Wigner
	RooRealVar *bwmean;
	RooRealVar *width;
	RooRealVar *bwmeangen;
	RooRealVar *widthgen;
	RooBreitWigner *bw;
	RooBreitWigner *bwgen;
	//This is CrystalBall
	RooRealVar *cbmean;
	RooRealVar *cbsigma;
	RooRealVar *cbalpha;
	RooRealVar *cbn;
	RooRealVar *cbmeangen;
	RooRealVar *cbsigmagen;
	RooRealVar *cbalphagen;
	RooRealVar *cbngen;
	RooCBShape *cb;
	RooCBShape *cbgen;
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
	RooFFTConvPdf* convpdfgen;

	double meanshift = 0;
	double meanshifterror = 0;

	double widthshift = 0;
	double widtherror = 0;

	Int_t h_normalized_data_entry[11] = {};

	//new roorealvar area:

	RooRealVar *roomass[11];
	Int_t cenlowlimit[11] = {0,10,20,30,50,0,15,70,0,14,0};
	Int_t cenhighlimit[11] = {10,20,30,50,100,15,100,90,14,100,100};






};

fitting_helper::fitting_helper(){
}

void fitting_helper::gethistograms(int opt,bool useunbin){

		if (opt == 1 || opt == 7)f1 = new TFile(raw_file, "READ");
		if (opt == 2)f1 = new TFile(ycut_file, "READ");
		if (opt == 3)f1 = new TFile(eff_file, "READ");
		if (opt == 4)f1 = new TFile(ycut_eff_file, "READ");
		if (opt == 5)f1 = new TFile(etacut_file,"READ");
		if (opt == 6)f1 = new TFile(etacut_eff_file,"READ");
		for (int i=0; i<11; i++){
			h_normalized_mc[i] = (TH1D*) f1->Get(Form("Normalized_mc_%i",i));
			h_normalized_data[i] = (TH1D*) f1->Get(Form("Normalized_data_%i",i));
			h_normalized_mc_bk[i] = (TH1D*) f1->Get(Form("Normalized_mc_bk_%i",i));
			h_normalized_data_entry[i] = h_normalized_data[i]->GetEntries();
		}

		f2 = new TFile(data_file,"READ");
		for (int i=0; i<11; i++){
			if (opt == 1) roodata[i] = (RooDataSet*)f2->Get(Form("roodata_%i",i));
			if (opt == 2) roodata[i] = (RooDataSet*)f2->Get(Form("roodata_y_%i",i));
			if (opt == 4) roodata[i] = (RooDataSet*)f2->Get(Form("roodata_y_eff_%i",i));
			if (opt == 3) roodata[i] = (RooDataSet*)f2->Get(Form("roodata_eff_%i",i));
			if (opt == 6) roodata[i] = (RooDataSet*)f2->Get(Form("roodata_eta_eff_%i",i));
			if (opt == 5) roodata[i] = (RooDataSet*)f2->Get(Form("roodata_eta_%i",i));
			if (opt == 7) roodata[i] = (RooDataSet*)f2->Get(Form("roodata_raw_uniform_%i",i));
		}
	//f2 = new TFile(mc_file_path,"READ");
}

void fitting_helper::setupfittingfunction(int cent, int type, bool useunbin){
	x = new RooRealVar("roomass","roomass",60,120);
	x->setBinning(RooBinning(10000,60,120),"cache");

	bwmean = new RooRealVar("bwmean", "bwmean",bwmeanip,80,100);
	width = new RooRealVar("width", "width",widthip,0,10 );
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width );

	/*if (type == 3 && useunbin){
		//raw + eff for unbinned
		cbmean = new RooRealVar("mean","mean",cbmeanip);
		cbsigma = new RooRealVar("sigma","sigma",sigmaip,0,10);
		cbalpha = new RooRealVar("alpha","alpha",1.7681,-5,5);
		cbn = new RooRealVar("n","n", 1.0567,-5,5);
		cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	}
	else if (type == 4 && useunbin){
		//y + eff for unbinned
		cbmean = new RooRealVar("mean","mean",cbmeanip);
		cbsigma = new RooRealVar("sigma","sigma",sigmaip,0,10);
		cbalpha = new RooRealVar("alpha","alpha",1.7934,-5,5);
		cbn = new RooRealVar("n","n", 1.0584,-5,5);
		cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	}
	else if (type == 6 && useunbin){
		//eta + eff for unbinned
		cbmean = new RooRealVar("mean","mean",cbmeanip);
		cbsigma = new RooRealVar("sigma","sigma",sigmaip,0,10);
		//cbsigma = new RooRealVar("sigma","sigma",0.00001);
		cbalpha = new RooRealVar("alpha","alpha",1.8529,-5,5);
		cbn = new RooRealVar("n","n", 0.8753,-5,5);
		cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	}
	else{*/
		cbmean = new RooRealVar("mean","mean",cbmeanip);
		//cbsigma = new RooRealVar("sigma","sigma",sigmaip,0,10);
		cbsigma = new RooRealVar("sigma","sigma",0.00001);
		cbalpha = new RooRealVar("alpha","alpha",alphaip,-5,5);
		cbn = new RooRealVar("n","n", nip ,-5,5);
		cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	//}

	decay = new RooRealVar("decay","decay",decayip,-5,+5);
	exp = new RooExponential("exp", "Exponential PDF", *x, *decay);

	newconvpdf = new RooFFTConvPdf("newconvpdf","newconvpdf", *x, *bw, *cb);

	fsig = new RooRealVar("fsig","signal fraction",fsigip,0,1);
    purepdf = new RooAddPdf("purepdf","cbconbw+exp",RooArgList(*newconvpdf,*exp),*fsig);


	background = new RooDataHist("background","background",*x,h_normalized_mc_bk[cent]);
	histpdf1 = new RooHistPdf("histpdf1", "histpdf1", *x, *background, 0);
	sig1frac = new RooRealVar("sig1frac", "fraction of component 1 in signal", sig1fracip,0,1);
	temp_sig = new RooAddPdf("temp_sig","temp_sig",RooArgList(*newconvpdf,*histpdf1),*sig1frac);

	data = new RooDataHist("data","data",*x,h_normalized_data[cent]);
	mc = new RooDataHist("mc","mc",*x,h_normalized_mc[cent]);


	/*for (int i = 0; i < 11; i++){
		roomass[i] = new RooRealVar(Form("roomass_%i",i),Form("roomass_%i",i),60,120);
		roomass[i]->setBinning(RooBinning(10000,60,120),"cache");


	}*/
	//pdf->setConvolutionWindow(*convocenter, *convowidth, 1);
	//RooFit::Minimizer("Minuit","migrad")
}

RooFitResult* fitting_helper::fitting(int opt1,int cent, bool useunbin,bool usew2unbin){
	//pdf->fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	//if (opt1 == 1 )pdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true));
	RooFitResult* fitresult;
	if (useunbin){
		if (!usew2unbin){
			if (opt1 == 1 )fitresult = newconvpdf->fitTo(*roodata[cent],RooFit::Save(true),RooFit::NumCPU(8),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
			if (opt1 == 2 )fitresult = purepdf->fitTo(*roodata[cent],RooFit::Save(true),RooFit::NumCPU(8),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
			if (opt1 == 3)fitresult = temp_sig->fitTo(*roodata[cent],RooFit::Save(true),RooFit::NumCPU(8),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
		}
		if (usew2unbin){
			if (opt1 == 1 )fitresult = newconvpdf->fitTo(*roodata[cent],RooFit::Save(true),RooFit::NumCPU(8),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (opt1 == 2 )fitresult = purepdf->fitTo(*roodata[cent],RooFit::Save(true),RooFit::NumCPU(8),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (opt1 == 3)fitresult = temp_sig->fitTo(*roodata[cent],RooFit::Save(true),RooFit::NumCPU(8),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
		}
	}

	if (!useunbin){
		if (opt1 == 1 )fitresult = newconvpdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
		if (opt1 == 2 )fitresult = purepdf->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
		if (opt1 == 3)fitresult = temp_sig->fitTo(*data,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
	}
	

	return fitresult;
}
void fitting_helper::plotting(int opt1, int opt, int cent, bool useunbin, bool usew2unbin){
	frame = x->frame(RooFit::Title("Z mass fit"));
	framecheck = x->frame(RooFit::Title("Background"));
	if (useunbin){
		if (opt == 1) {
		roodata[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(60),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		}
		if (opt == 2) {
		roodata[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(60),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		purepdf->plotOn(frame,RooFit::Components(*exp), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		purepdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		purepdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		purepdf->plotOn(framecheck,RooFit::Components(*exp), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		}
		if (opt == 3) {
		roodata[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(60),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		temp_sig->plotOn(frame,RooFit::Components(*histpdf1), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		temp_sig->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		temp_sig->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		temp_sig->plotOn(framecheck,RooFit::Components(*histpdf1), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		}
	}
	if (!useunbin){
		if (opt == 1) {
		data->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(60),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		}
		if (opt == 2) {
		data->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(60),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		purepdf->plotOn(frame,RooFit::Components(*exp), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		purepdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		purepdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		purepdf->plotOn(framecheck,RooFit::Components(*exp), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		}
		if (opt == 3) {
		data->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(60),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		temp_sig->plotOn(frame,RooFit::Components(*histpdf1), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		temp_sig->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		temp_sig->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		temp_sig->plotOn(framecheck,RooFit::Components(*histpdf1), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		}
	}
	

	//Create pull plot:
	if (useunbin){
		if (opt ==1){
			residuals = frame->residHist("roodata","fit",true,false); //true = pull, false = center of the bin
		}
		if (opt ==2)residuals = frame->residHist("roodata","fit",true,false); //true = pull, false = center of the bin
		if (opt ==3)residuals = frame->residHist("roodata","fit",true,false); //true = pull, false = center of the bin

	}
	if (!useunbin){
		if (opt ==1)residuals = frame->residHist("roodata","fit",true,false); //true = pull, false = center of the bin
		if (opt ==2)residuals = frame->residHist("roodata","fit",true,false); //true = pull, false = center of the bin
		if (opt ==3)residuals = frame->residHist("roodata","fit",true,false); //true = pull, false = center of the bin

	}
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
	if (useunbin && usew2unbin){
		if (opt1 == 1) textBox->AddText("Raw,w2");
		if (opt1 == 2) textBox->AddText("Zy<1,w2");
		if (opt1 == 3) textBox->AddText("Eff,w2");
		if (opt1 == 4) textBox->AddText("Zy<1+Eff,w2");
		if (opt1 == 5) textBox->AddText("ueta<1,w2");
		if (opt1 == 6) textBox->AddText("ueta<1+Eff,w2");
		if (opt1 == 7) textBox->AddText("Raw_uniform,w2");
	}
	if (useunbin && !usew2unbin){
		if (opt1 == 1) textBox->AddText("Raw,Asym");
		if (opt1 == 2) textBox->AddText("Zy<1,Asym");
		if (opt1 == 3) textBox->AddText("Eff,Asym");
		if (opt1 == 4) textBox->AddText("Zy<1+Eff,Asym");
		if (opt1 == 5) textBox->AddText("ueta<1,Asym");
		if (opt1 == 6) textBox->AddText("ueta<1+Eff,Asym");
		if (opt1 == 7) textBox->AddText("Raw_uniform,Asym");
	}
	if (!useunbin){
		if (opt1 == 1) textBox->AddText("Raw,binned");
		if (opt1 == 2) textBox->AddText("Zy<1,binned");
		if (opt1 == 3) textBox->AddText("Eff,binned");
		if (opt1 == 4) textBox->AddText("Zy<1+Eff,binned");
		if (opt1 == 5) textBox->AddText("ueta<1,binned");
		if (opt1 == 6) textBox->AddText("ueta<1+Eff,binned");
		if (opt1 == 7) textBox->AddText("Raw_uniform,binned");
	}
    textBox->Draw();
	c1->cd(2);
	pullFrame->Draw();
	c1->cd(3);
	framecheck->Draw();
	if (useunbin){
		if (!usew2unbin){
			if (opt1 == 1){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/unbin/raw/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/unbin/raw/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/unbin/raw/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 2){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/unbin/ycut/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/unbin/ycut/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/unbin/ycut/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 3){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/unbin/eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/unbin/eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/unbin/eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 4){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/unbin/ycut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/unbin/ycut_eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/unbin/ycut_eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 5){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/unbin/etacut/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/unbin/etacut/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/unbin/etacut/cb+bw+temp_%i.pdf",cent));

	}
	if (opt1 == 6){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/unbin/etacut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/unbin/etacut_eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/unbin/etacut_eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 7){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/unbin/raw_uniform/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/unbin/raw_uniform/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/unbin/raw_uniform/cb+bw+temp_%i.pdf",cent));
	}
	}

	if (usew2unbin){
		if (opt1 == 1){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/w2unbin/raw/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/w2unbin/raw/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/w2unbin/raw/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 2){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/w2unbin/ycut/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/w2unbin/ycut/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/w2unbin/ycut/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 3){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/w2unbin/eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/w2unbin/eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/w2unbin/eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 4){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/w2unbin/ycut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/w2unbin/ycut_eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/w2unbin/ycut_eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 5){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/w2unbin/etacut/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/w2unbin/etacut/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/w2unbin/etacut/cb+bw+temp_%i.pdf",cent));

	}
	if (opt1 == 6){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/w2unbin/etacut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/w2unbin/etacut_eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/w2unbin/etacut_eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 7){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/w2unbin/raw_uniform/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/w2unbin/raw_uniform/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/w2unbin/raw_uniform/cb+bw+temp_%i.pdf",cent));
	}
	}
	}

	if (!useunbin){
		if (opt1 == 1){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/binned/raw/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/binned/raw/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/binned/raw/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 2){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/binned/ycut/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/binned/ycut/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/binned/ycut/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 3){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/binned/eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/binned/eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/binned/eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 4){
		if (opt == 1) c1->SaveAs(Form("./fitresultplot/binned/ycut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2) c1->SaveAs(Form("./fitresultplot/binned/ycut_eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3) c1->SaveAs(Form("./fitresultplot/binned/ycut_eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 5){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/binned/etacut/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/binned/etacut/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/binned/etacut/cb+bw+temp_%i.pdf",cent));

	}
	if (opt1 == 6){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/binned/etacut_eff/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/binned/etacut_eff/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/binned/etacut_eff/cb+bw+temp_%i.pdf",cent));
	}
	if (opt1 == 7){
		if (opt == 1)c1->SaveAs(Form("./fitresultplot/binned/raw_uniform/cb+bw_%i.pdf",cent));
		if (opt == 2)c1->SaveAs(Form("./fitresultplot/binned/raw_uniform/cb+bw+exp_%i.pdf",cent));
		if (opt == 3)c1->SaveAs(Form("./fitresultplot/binned/raw_uniform/cb+bw+temp_%i.pdf",cent));
	}
	}


}

void fitting_helper::setupfittingfunctionchi(){
	meanshift = 2;
}

void fitting_helper::plotTGraph(TGraphErrors* x1, TGraphErrors* x2, TGraphErrors *x3, int opt, int opt1,bool useunbin,bool usew2unbin){
	gStyle->SetEndErrorSize(4);
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
	
	if (opt == 1) x1->GetYaxis()->SetRangeUser(-0.6,0.2);
	if (opt == 2) x1->GetYaxis()->SetRangeUser(-1,1.4);

	x1->SetLineWidth(2);
	x2->SetLineWidth(2);
	x3->SetLineWidth(2);

	x1->SetMarkerSize(1.3);
	x1->SetMarkerStyle(25);
	x2->SetMarkerSize(1.3);
	x2->SetMarkerStyle(4);
	x3->SetMarkerSize(1.3);
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

	TLegend *legend = new TLegend(0.6, 0.65, 0.8, 0.90);
	legend->SetTextFont(40);
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
	if (useunbin && usew2unbin){
		if (opt1 == 1) pt->AddText("Raw,w2");
		if (opt1 == 2) pt->AddText("Zy<1,w2");
		if (opt1 == 3) pt->AddText("Eff,w2");
		if (opt1 == 4) pt->AddText("Zy<1+Eff,w2");
		if (opt1 == 5) pt->AddText("ueta<1,w2");
		if (opt1 == 6) pt->AddText("ueta<1+Eff,w2");
		if (opt1 == 7) pt->AddText("Raw_uniform,w2");
	}
	if (useunbin && !usew2unbin){
		if (opt1 == 1) pt->AddText("Raw,Asym");
		if (opt1 == 2) pt->AddText("Zy<1,Asym");
		if (opt1 == 3) pt->AddText("Eff,Asym");
		if (opt1 == 4) pt->AddText("Zy<1+Eff,Asym");
		if (opt1 == 5) pt->AddText("ueta<1,Asym");
		if (opt1 == 6) pt->AddText("ueta<1+Eff,Asym");
		if (opt1 == 7) pt->AddText("Raw_uniform,Asym");
	}
	if (!useunbin){
		if (opt1 == 1) pt->AddText("Raw,binned");
		if (opt1 == 2) pt->AddText("Zy<1,binned");
		if (opt1 == 3) pt->AddText("Eff,binned");
		if (opt1 == 4) pt->AddText("Zy<1+Eff,binned");
		if (opt1 == 5) pt->AddText("ueta<1,binned");
		if (opt1 == 6) pt->AddText("ueta<1+Eff,binned");
		if (opt1 == 7) pt->AddText("Raw_uniform,binned");
	}

	legend->Draw();
	pt->Draw();

	if (useunbin){
		if (!usew2unbin){
			if (opt == 1 && opt1 == 1)	c_2->SaveAs("./dM_dSig/unbin/raw/dM.pdf");
			if (opt == 2 && opt1 == 1) 	c_2->SaveAs("./dM_dSig/unbin/raw/dSig.pdf");

			if (opt == 1 && opt1 == 2)	c_2->SaveAs("./dM_dSig/unbin/ycut/dM.pdf");
			if (opt == 2 && opt1 == 2) 	c_2->SaveAs("./dM_dSig/unbin/ycut/dSig.pdf");

			if (opt == 1 && opt1 == 3)	c_2->SaveAs("./dM_dSig/unbin/eff/dM.pdf");
			if (opt == 2 && opt1 == 3) 	c_2->SaveAs("./dM_dSig/unbin/eff/dSig.pdf");

			if (opt == 1 && opt1 == 4)	c_2->SaveAs("./dM_dSig/unbin/ycut_eff/dM.pdf");
			if (opt == 2 && opt1 == 4) 	c_2->SaveAs("./dM_dSig/unbin/ycut_eff/dSig.pdf");

			if (opt == 1 && opt1 == 5)	c_2->SaveAs("./dM_dSig/unbin/etacut/dM.pdf");
			if (opt == 2 && opt1 == 5) 	c_2->SaveAs("./dM_dSig/unbin/etacut/dSig.pdf");

			if (opt == 1 && opt1 == 6)	c_2->SaveAs("./dM_dSig/unbin/etacut_eff/dM.pdf");
			if (opt == 2 && opt1 == 6) 	c_2->SaveAs("./dM_dSig/unbin/etacut_eff/dSig.pdf");

			if (opt == 1 && opt1 == 7)	c_2->SaveAs("./dM_dSig/unbin/raw_uniform/dM.pdf");
			if (opt == 2 && opt1 == 7) 	c_2->SaveAs("./dM_dSig/unbin/raw_uniform/dSig.pdf");
		}
		if (usew2unbin){
			if (opt == 1 && opt1 == 1)	c_2->SaveAs("./dM_dSig/w2unbin/raw/dM.pdf");
			if (opt == 2 && opt1 == 1) 	c_2->SaveAs("./dM_dSig/w2unbin/raw/dSig.pdf");

			if (opt == 1 && opt1 == 2)	c_2->SaveAs("./dM_dSig/w2unbin/ycut/dM.pdf");
			if (opt == 2 && opt1 == 2) 	c_2->SaveAs("./dM_dSig/w2unbin/ycut/dSig.pdf");

			if (opt == 1 && opt1 == 3)	c_2->SaveAs("./dM_dSig/w2unbin/eff/dM.pdf");
			if (opt == 2 && opt1 == 3) 	c_2->SaveAs("./dM_dSig/w2unbin/eff/dSig.pdf");

			if (opt == 1 && opt1 == 4)	c_2->SaveAs("./dM_dSig/w2unbin/ycut_eff/dM.pdf");
			if (opt == 2 && opt1 == 4) 	c_2->SaveAs("./dM_dSig/w2unbin/ycut_eff/dSig.pdf");

			if (opt == 1 && opt1 == 5)	c_2->SaveAs("./dM_dSig/w2unbin/etacut/dM.pdf");
			if (opt == 2 && opt1 == 5) 	c_2->SaveAs("./dM_dSig/w2unbin/etacut/dSig.pdf");

			if (opt == 1 && opt1 == 6)	c_2->SaveAs("./dM_dSig/w2unbin/etacut_eff/dM.pdf");
			if (opt == 2 && opt1 == 6) 	c_2->SaveAs("./dM_dSig/w2unbin/etacut_eff/dSig.pdf");

			if (opt == 1 && opt1 == 7)	c_2->SaveAs("./dM_dSig/w2unbin/raw_uniform/dM.pdf");
			if (opt == 2 && opt1 == 7) 	c_2->SaveAs("./dM_dSig/w2unbin/raw_uniform/dSig.pdf");

		}
	}

	if (!useunbin){
	if (opt == 1 && opt1 == 1)	c_2->SaveAs("./dM_dSig/binned/raw/dM.pdf");
	if (opt == 2 && opt1 == 1) 	c_2->SaveAs("./dM_dSig/binned/raw/dSig.pdf");

	if (opt == 1 && opt1 == 2)	c_2->SaveAs("./dM_dSig/binned/ycut/dM.pdf");
	if (opt == 2 && opt1 == 2) 	c_2->SaveAs("./dM_dSig/binned/ycut/dSig.pdf");

	if (opt == 1 && opt1 == 3)	c_2->SaveAs("./dM_dSig/binned/eff/dM.pdf");
	if (opt == 2 && opt1 == 3) 	c_2->SaveAs("./dM_dSig/binned/eff/dSig.pdf");

	if (opt == 1 && opt1 == 4)	c_2->SaveAs("./dM_dSig/binned/ycut_eff/dM.pdf");
	if (opt == 2 && opt1 == 4) 	c_2->SaveAs("./dM_dSig/binned/ycut_eff/dSig.pdf");

	if (opt == 1 && opt1 == 5)	c_2->SaveAs("./dM_dSig/binned/etacut/dM.pdf");
	if (opt == 2 && opt1 == 5) 	c_2->SaveAs("./dM_dSig/binned/etacut/dSig.pdf");

	if (opt == 1 && opt1 == 6)	c_2->SaveAs("./dM_dSig/binned/etacut_eff/dM.pdf");
	if (opt == 2 && opt1 == 6) 	c_2->SaveAs("./dM_dSig/binned/etacut_eff/dSig.pdf");

	if (opt == 1 && opt1 == 7)	c_2->SaveAs("./dM_dSig/binned/raw_uniform/dM.pdf");
	if (opt == 2 && opt1 == 7) 	c_2->SaveAs("./dM_dSig/binned/raw_uniform/dSig.pdf");	
	}



}

void fitting_helper::getmcdataset(){
	f1 = new TFile(mc_file,"READ");
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
}

void fitting_helper::setupfittingfunctionmc(){

	x = new RooRealVar("roomass","roomass",60,120);
	x->setBinning(RooBinning(10000,60,120),"cache");

	bwmean = new RooRealVar("bwmean", "bwmean",bwmeanrecoip,80,100);
	width = new RooRealVar("width", "width",widthrecoip,0,10 );
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width );
	bwmeangen = new RooRealVar("bwmeangen","bwmeangen",bwmeangenip,80,100);
	widthgen = new RooRealVar("widthgen","widthgen",widthgenip,0,10);
	bwgen = new RooBreitWigner("bwgen","bwgen",*x,*bwmeangen,*widthgen);

	cbmean = new RooRealVar("mean","mean",cbmeanrecoip);
	//cbsigma = new RooRealVar("sigma","sigma",sigmarecoip,0,10); // Here Fixed sig for cb, subject to change
	cbsigma = new RooRealVar("sigma","sigma",0.00001);
	cbalpha = new RooRealVar("alpha","alpha",alpharecoip,-5,5);
	cbn = new RooRealVar("n","n", nrecoip ,-5,5);
	cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);

	cbmeangen = new RooRealVar("cbmeangen","cbmeangen",cbmeangenip);
	cbsigmagen = new RooRealVar("cbsigmagen","cbsigmagen",sigmagenip,0,10);
	cbalphagen = new RooRealVar("cbalphagen","cbalphagen",alphagenip,-5,5);
	cbngen = new RooRealVar("cbngen","cbngen",ngenip,-5,5);
	cbgen = new RooCBShape("cbgen","cbgen",*x,*cbmeangen,*cbsigmagen,*cbalphagen,*cbngen);


	newconvpdf = new RooFFTConvPdf("newconvpdf","newconvpdf", *x, *bw, *cb); // This is CB + BW
	convpdfgen = new RooFFTConvPdf("convpdfgen","convpdfgen",*x,*bwgen,*cbgen);

}

RooFitResult* fitting_helper::fittingmc(int type, bool iseff, bool usew2unbin, int cent){
	RooFitResult* fitresult;
	if (type == 1){
		if (iseff){
			if (usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_eff[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (!usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_eff[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
		}
		if (!iseff){
			if (usew2unbin)fitresult = newconvpdf->fitTo(*rooreco[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (!usew2unbin)fitresult = newconvpdf->fitTo(*rooreco[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
		}
	}
	if (type == 2){
		if (iseff){
			if (usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_y_eff[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (!usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_y_eff[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
		}
		if (!iseff){
			if (usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_y[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (!usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_y[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
		}
	}
	if (type == 3){
		if (iseff){
			if (usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_eta_eff[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (!usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_eta_eff[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
		}
		if (!iseff){
			if (usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_eta[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
			if (!usew2unbin)fitresult = newconvpdf->fitTo(*rooreco_eta[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
		}
	}
	if (type == 4){
		if (usew2unbin)fitresult = convpdfgen->fitTo(*roogen[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
		if (!usew2unbin)fitresult = convpdfgen->fitTo(*roogen[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
	}
	if (type == 5){
		if (usew2unbin)fitresult = convpdfgen->fitTo(*roogen_y[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
		if (!usew2unbin)fitresult = convpdfgen->fitTo(*roogen_y[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));
	}
	if (type == 6){
		if (usew2unbin)fitresult = convpdfgen->fitTo(*roogen_eta[cent],RooFit::Save(true),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
		if (!usew2unbin)fitresult = convpdfgen->fitTo(*roogen_eta[cent],RooFit::Save(true),RooFit::AsymptoticError(true),RooFit::Minimizer("Minuit2","migrad"));	
	}
	return fitresult;
}

void fitting_helper::plottingmc(int type, int cent, bool iseff, bool usew2unbin){
	frame = x->frame(RooFit::Title("Z mass fit"));
	Int_t numEntries;
	Double_t integral;
	if (type == 1){
		if (iseff){
			rooreco_eff[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
			newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
			numEntries = rooreco_eff[cent]->numEntries();
			integral = rooreco_eff[cent]->sumEntries();
		}
		if (!iseff){
			rooreco[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
			newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
			numEntries = rooreco[cent]->numEntries();
			integral = rooreco[cent]->sumEntries();
		}
	}
	if (type == 2){
		if (iseff){
			rooreco_y_eff[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
			newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
			numEntries = rooreco_y_eff[cent]->numEntries();
			integral = rooreco_y_eff[cent]->sumEntries();
		}
		if (!iseff){
			rooreco_y[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
			newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
			numEntries = rooreco_y[cent]->numEntries();
			integral = rooreco_y[cent]->sumEntries();
		}
	}
	if (type == 3){
		if (iseff){
			rooreco_eta_eff[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
			newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
			numEntries = rooreco_eta_eff[cent]->numEntries();
			integral = rooreco_eta_eff[cent]->sumEntries();
		}
		if (!iseff){
			rooreco_eta[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
			newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
			numEntries = rooreco_eta[cent]->numEntries();
			integral = rooreco_eta[cent]->sumEntries();
		}
	}
	if (type == 4){
		roogen[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		convpdfgen->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		convpdfgen->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
		numEntries = roogen[cent]->numEntries();
		integral = roogen[cent]->sumEntries();
	}
	if (type == 5){
		roogen_y[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		convpdfgen->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		convpdfgen->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
		numEntries = roogen_y[cent]->numEntries();
		integral = roogen_y[cent]->sumEntries();
		
	}
	if (type == 6){
		roogen_eta[cent]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		convpdfgen->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
		convpdfgen->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.5,1,0.9),RooFit::ShowConstants(kTRUE));
		numEntries = roogen_eta[cent]->numEntries();
		integral = roogen_eta[cent]->sumEntries();
	}

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

	double meanofBW;
	double meanerrorofBW;
	double widthofBW;
	double widtherrorofBW;

	if (type == 1 || type == 2 || type == 3){
		meanofBW = bwmean->getVal() - 91.1876;
		meanerrorofBW = bwmean->getError();
		widthofBW = width->getVal() - 2.4955;
		widtherrorofBW = width->getError();

		meanshift = meanofBW;
		meanshifterror = meanerrorofBW;
		widthshift = widthofBW;
		widtherror = widtherrorofBW;
	}
	if (type == 4 || type == 5 || type == 6){
		meanofBW = bwmeangen->getVal() - 91.1876;
		meanerrorofBW = bwmeangen->getError();
		widthofBW = widthgen->getVal() - 2.4955;
		widtherrorofBW = widthgen->getError();

		meanshift = meanofBW;
		meanshifterror = meanerrorofBW;
		widthshift = widthofBW;
		widtherror = widtherrorofBW;
	}


	textBox = new TPaveText(0.1, 0.5, 0.4, 0.9, "NDC");
    textBox->SetFillColor(0);
    textBox->SetTextSize(0.03);
    textBox->AddText(Form("Mean value: %.2f", meanofBW));
	textBox->AddText(Form("Width value: %.2f",widthofBW));
	textBox->AddText(Form("# of Entries: %i",numEntries));
	textBox->AddText(Form("Integral: %.2f",integral));
	if (usew2unbin){
		if (type == 1){
			if (iseff) textBox->AddText("Reco_raw_eff,w2");
			if (!iseff) textBox->AddText("Reco_raw,w2");
		}
		if (type == 2){
			if (iseff) textBox->AddText("reco_y_eff,w2");
			if (!iseff) textBox->AddText("reco_y,w2");
		}
		if (type == 3){
			if (iseff) textBox->AddText("reco_eta_eff,w2");
			if (!iseff) textBox->AddText("reco_eta,w2");
		}
		if (type == 4){
			textBox->AddText("gen_raw,w2");
		}
		if (type == 5){
			textBox->AddText("gen_y,w2");
		}
		if (type == 6){
			textBox->AddText("gen_eta,w2");
		}
	}
	if (!usew2unbin){
		if (type == 1){
			if (iseff) textBox->AddText("Reco_raw_eff,Asym");
			if (!iseff) textBox->AddText("Reco_raw,Asym");
		}
		if (type == 2){
			if (iseff) textBox->AddText("reco_y_eff,Asym");
			if (!iseff) textBox->AddText("reco_y,Asym");
		}
		if (type == 3){
			if (iseff) textBox->AddText("reco_eta_eff,Asym");
			if (!iseff) textBox->AddText("reco_eta,Asym");
		}
		if (type == 4){
			textBox->AddText("gen_raw,Asym");
		}
		if (type == 5){
			textBox->AddText("gen_y,Asym");
		}
		if (type == 6){
			textBox->AddText("gen_eta,Asym");
		}

	}
    textBox->Draw();	

	c1->cd(2);
	pullFrame->Draw();

	if (type == 1){
		if (iseff){
			if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/reco_raw_eff/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/reco_raw_eff/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		}
		if (!iseff){
			if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/reco_raw/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/reco_raw/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		}
	}
	if (type == 2){
		if (iseff){
			if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/reco_y_eff/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/reco_y_eff/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		}
		if (!iseff){
			if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/reco_y/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/reco_y/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		}
	}
	if (type == 3){
		if (iseff){
			if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/reco_eta_eff/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/reco_eta_eff/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		}
		if (!iseff){
			if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/reco_eta/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/reco_eta/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		}
	}
	if (type == 4){
		if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/gen_raw/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/gen_raw/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
	}
	if (type == 5){
		if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/gen_y/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/gen_y/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		
	}
	if (type == 6){
		if (usew2unbin)c1->SaveAs(Form("./fitresultmc/w2unbin/gen_eta/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		if (!usew2unbin)c1->SaveAs(Form("./fitresultmc/gen_eta/cb+bw_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
	}



}

void fitting_helper::plotTGraphmc(TGraphErrors* x1, int type, bool iseff, bool usew2unbin, int param, int cent){
	gStyle->SetEndErrorSize(4);
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

	double cbsigmaplot = cbsigma->getVal();
	double cbsigmaerror = cbsigma->getError();
	double cbalphaplot = cbalpha->getVal();
	double cbalphaerror = cbalpha->getError();
	double cbnplot = cbn->getVal();
	double cbnerror = cbn->getError();

	TCanvas *c2 = new TCanvas("","",600,600);
	TPad *pad = new TPad("pad", "Pad", 0, 0, 1, 1);
    pad->Draw();
    pad->cd();
	x1->SetTitle("");
	x1->GetXaxis()->SetTitle("Central centrality");
	x1->GetXaxis()->CenterTitle();
	if (param == 1) x1->GetYaxis()->SetTitle("dM(GeV)");
	if (param == 2) x1->GetYaxis()->SetTitle("dWidth(GeV)");
	if (param == 3) x1->GetYaxis()->SetTitle("Sigma");
	if (param == 4)	x1->GetYaxis()->SetTitle("alpha");
	if (param == 5)	x1->GetYaxis()->SetTitle("n");
	x1->GetYaxis()->CenterTitle();
	x1->GetYaxis()->SetTitleOffset(1.9);

	if (param == 1) x1->GetYaxis()->SetRangeUser(-0.2,0.1);
	if (param == 2) x1->GetYaxis()->SetRangeUser(-0.2,0.4);
	if (param == 3) x1->GetYaxis()->SetRangeUser(0.6,1.3);
	if (param == 4)	x1->GetYaxis()->SetRangeUser(1.6,2);
	if (param == 5)	x1->GetYaxis()->SetRangeUser(0.8,1.3);

	x1->SetLineWidth(2);
	x1->SetMarkerSize(1.5);
	x1->SetMarkerStyle(25);
	x1->SetMarkerColor(2);
	TLegend *legend = new TLegend(0.6, 0.7, 0.8, 0.9);
	legend->SetTextFont(40);
	legend->SetBorderSize(0);
    legend->AddEntry(x1, "cb+bw", "lep"); 
	legend->SetFillColor(kWhite);
	x1->Draw("AP");
	//legend->Draw();

	TF1 *f1 = new TF1("f1", "[0]", 0, 80);

	f1->SetLineColor(2);

	x1->Fit(f1,"QR");

	double fittedValue = f1->GetParameter(0);
	double error = f1->GetParError(0);

	TPaveText *pt = new TPaveText(0.2, 0.8, 0.55, 0.9, "brNDC"); // normalized coordinates
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextAlign(12); // Align left and vertically centered
	pt->SetTextFont(42);
	pt->SetTextSize(0.02);
	pt->SetMargin(0.02);
	pt->AddText(Form("Fitted constant for cb+bw = %.4f #pm %.4f", fittedValue, error));
	
	if (usew2unbin){
		if (type == 1){
			if (iseff) pt->AddText("Reco_raw_eff,w2");
			if (!iseff) pt->AddText("Reco_raw,w2");
		}
		if (type == 2){
			if (iseff) pt->AddText("reco_y_eff,w2");
			if (!iseff) pt->AddText("reco_y,w2");
		}
		if (type == 3){
			if (iseff) pt->AddText("reco_eta_eff,w2");
			if (!iseff) pt->AddText("reco_eta,w2");
		}
		if (type == 4){
			pt->AddText("gen_raw,w2");
		}
		if (type == 5){
			pt->AddText("gen_y,w2");
		}
		if (type == 6){
			pt->AddText("gen_eta,w2");
		}
	}
	if (!usew2unbin){
		if (type == 1){
			if (iseff) pt->AddText("Reco_raw_eff,Asym");
			if (!iseff) pt->AddText("Reco_raw,Asym");
		}
		if (type == 2){
			if (iseff) pt->AddText("reco_y_eff,Asym");
			if (!iseff) pt->AddText("reco_y,Asym");
		}
		if (type == 3){
			if (iseff) pt->AddText("reco_eta_eff,Asym");
			if (!iseff) pt->AddText("reco_eta,Asym");
		}
		if (type == 4){
			pt->AddText("gen_raw,Asym");
		}
		if (type == 5){
			pt->AddText("gen_y,Asym");
		}
		if (type == 6){
			pt->AddText("gen_eta,Asym");
		}
	}

	pt->Draw();

	if (usew2unbin){
		if (type == 1){
			if (iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw_eff/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw_eff/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw_eff/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw_eff/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw_eff/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
			if (!iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_raw/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
		}
		if (type == 2){
			if (iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y_eff/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y_eff/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y_eff/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y_eff/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y_eff/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
			if (!iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_y/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}

		}
		if (type == 3){
			if (iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta_eff/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta_eff/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta_eff/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta_eff/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta_eff/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
			if (!iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/reco_eta/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}

		}
		if (type == 4){
			if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_raw/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_raw/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_raw/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_raw/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_raw/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
		}
		if (type == 5){
			if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_y/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_y/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_y/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_y/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_y/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));

		}
		if (type == 6){
			if (param == 1)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_eta/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 2)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_eta/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 3)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_eta/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 4)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_eta/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 5)c2->SaveAs(Form("./parameters_mc/w2unbin/gen_eta/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));

		}
	}
	if (!usew2unbin){
		if (type == 1){
			if (iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw_eff/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw_eff/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw_eff/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw_eff/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw_eff/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
			if (!iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_raw/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
		}
		if (type == 2){
			if (iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y_eff/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y_eff/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y_eff/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y_eff/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y_eff/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
			if (!iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_y/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}

		}
		if (type == 3){
			if (iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta_eff/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta_eff/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta_eff/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta_eff/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta_eff/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}
			if (!iseff){
				if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
				if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/reco_eta/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			}

		}
		if (type == 4){
			if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_raw/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_raw/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_raw/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_raw/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_raw/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));

		}
		if (type == 5){
			if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_y/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_y/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_y/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_y/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_y/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));

		}
		if (type == 6){
			if (param == 1)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_eta/dM_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 2)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_eta/dWidth_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 3)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_eta/Sigma_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 4)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_eta/alpha_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));
			if (param == 5)c2->SaveAs(Form("./parameters_mc/asymunbin/gen_eta/N_%i_%i.pdf",cenlowlimit[cent],cenhighlimit[cent]));

		}

	}

}