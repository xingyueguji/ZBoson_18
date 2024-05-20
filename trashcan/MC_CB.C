#include "MC_18.h"

void MC_CB(int opt = 2, int opt1 = 2, int opt2 = 3){
    // opt = 1 means CB+BW opt = 2 means CB+BW+EXP opt = 3 means CB+BW+TEMPLATE
	// opt1 = 1 means MC, opt1 = 2 means data
	// opt2 == 1 means all, opt2 == 2 means 24, opt2 == 3 means 100

	TFile *f1 = new TFile("./mc_bk_plot/W_histogram.root","READ");
	TH1D* h_1;
	TH1D* h_2;
	TH1D* h_3;
	TH1D* h_4;
	TH1D* h_5;
	TH1D* h_6;
	TH1D* h_7;
	TH1D* h_8;
	TH1D* h_9;
	TH1D* h_10;

	if (opt2 == 1){

	h_1 = (TH1D*) f1->Get("ZGEN");
	h_2 = (TH1D*) f1->Get("ZTOUU");
	h_3 = (TH1D*) f1->Get("Tau_BK");
	h_4 = (TH1D*) f1->Get("mc_W_bk");
	h_5 = (TH1D*) f1->Get("mc_bk_tt");
	h_6 = (TH1D*) f1->Get("data_same_hist");
	h_7 = (TH1D*) f1->Get("data_hist");

	h_8 = (TH1D*) f1->Get("Normalized_data");
	h_9 = (TH1D*) f1->Get("Normalized_mc");
	h_10 = (TH1D*) f1->Get("Normalized_mc_bk");

	}
	if (opt2 == 2){

	h_1 = (TH1D*) f1->Get("ZGEN");
	h_2 = (TH1D*) f1->Get("ZTOUU_24");
	h_3 = (TH1D*) f1->Get("Tau_BK_24");
	h_4 = (TH1D*) f1->Get("mc_W_bk_24");
	h_5 = (TH1D*) f1->Get("mc_bk_tt_24");
	h_6 = (TH1D*) f1->Get("data_same_hist_24");
	h_7 = (TH1D*) f1->Get("data_hist_24");

	h_8 = (TH1D*) f1->Get("Normalized_data_24");
	h_9 = (TH1D*) f1->Get("Normalized_mc_24");
	h_10 = (TH1D*) f1->Get("Normalized_mc_bk_24");
		
	}
	if (opt2 == 3){

	h_1 = (TH1D*) f1->Get("ZGEN");
	h_2 = (TH1D*) f1->Get("ZTOUU_100");
	h_3 = (TH1D*) f1->Get("Tau_BK_100");
	h_4 = (TH1D*) f1->Get("mc_W_bk_100");
	h_5 = (TH1D*) f1->Get("mc_bk_tt_100");
	h_6 = (TH1D*) f1->Get("data_same_hist_100");
	h_7 = (TH1D*) f1->Get("data_hist_100");

	h_8 = (TH1D*) f1->Get("Normalized_data_100");
	h_9 = (TH1D*) f1->Get("Normalized_mc_100");
	h_10 = (TH1D*) f1->Get("Normalized_mc_bk_100");
		
	}
	/*cout << "Underflow is " << endl;
	cout << h_1->GetBinContent(121) << endl; // this is non-zero but okay because I did not use gen and the cut messed up with gen calculation (not 60 to 120)
	cout << h_2->GetBinContent(121) << endl;
	cout << h_3->GetBinContent(121) << endl;
	cout << h_4->GetBinContent(121) << endl;
	cout << h_5->GetBinContent(121) << endl;
	cout << h_6->GetBinContent(121) << endl;
	cout << h_7->GetBinContent(121) << endl;
	cout << h_8->GetBinContent(121) << endl;
	cout << h_9->GetBinContent(121) << endl;
	cout << h_10->GetBinContent(121) << endl;*/


	double bwmeanip = 90.892, widthip = 3.43133, cbmeanip = 0, sigmaip = 0.544, alphaip = 1.828, nip = 0.9378, decayip = -0.07613, fsigip = 0.99, sig1fracip = 0.99;

	//if(opt1 == 2) bwmeanip = 91, widthip = 4, cbmeanip = 0, sigmaip = 5.03, alphaip = 2, nip = 1, decayip = 0, fsigip = 0.99, sig1fracip = 0.99;

	//This is Breit-Wigner
	RooRealVar x("x", "x",60., 120. );
	RooRealVar bwmean("bwmean", "bwmean",bwmeanip,80,100);
	//RooRealVar m0("m0", "m0",m0ip,80,100);
	RooRealVar width("width", "width",widthip,0,10 );
	RooBreitWigner bw("bw", "bw", x, bwmean, width );

    //This is CrystalBall
	RooRealVar cbmean("mean","mean",cbmeanip);
	RooRealVar cbsigma("sigma","sigma",sigmaip,0,10);
	RooRealVar cbalpha("alpha","alpha",alphaip,-5,5);
	RooRealVar cbn("n","n", nip ,-5,5);
	RooCBShape cb("cb","cb", x, cbmean, cbsigma, cbalpha, cbn);
	//RooCBShape cb("cb","cb", x, m0, sigma, alpha, n );
    
	//This is Exponential

		RooRealVar decay("decay","decay",decayip,-5,5);
		RooExponential exp("exp", "Exponential PDF", x, decay);


    
	//This is BW + CB
	RooNumConvPdf pdf("pdf","pdf", x, bw, cb);

    //This is BW+CB+EXP

		RooRealVar fsig("fsig","signal fraction",fsigip,0,1);
		RooAddPdf purepdf("purepdf","cbconbw+exp",RooArgList(pdf,exp),fsig);


	//This is BW+CB+TEMPLATE

		RooDataHist background("background","background",x, h_10);
		RooDataHist *ptr_background = &background;
		RooHistPdf histpdf1("histpdf1", "histpdf1", x, *ptr_background, 0);
		RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", sig1fracip,0,1);
		RooAddPdf temp_sig("temp_sig","temp_sig",RooArgList(pdf,histpdf1),sig1frac);

	//RooDataHist data("data", "data", x, h_2);
	RooDataHist data("data","data",x, h_8);
	//RooDataHist mc("mc","mc",x,h_2);
	RooDataHist mc("mc","mc",x,h_9);	

	RooFitResult* result;

	if (opt == 1 && opt1 == 1) result = pdf.fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true),RooFit::MaxCalls(5000));
	if (opt == 1 && opt1 == 2) result = pdf.fitTo(data,RooFit::Save(true),RooFit::SumW2Error(true));
	if (opt == 2 && opt1 == 1) result = purepdf.fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true));
	if (opt == 2 && opt1 == 2) result = purepdf.fitTo(data,RooFit::Save(true),RooFit::SumW2Error(true));
	if (opt == 3 && opt1 == 1) result = temp_sig.fitTo(mc,RooFit::Save(true),RooFit::SumW2Error(true));
	if (opt == 3 && opt1 == 2) result = temp_sig.fitTo(data,RooFit::Save(true),RooFit::SumW2Error(true));

	
	
	result->floatParsFinal().Print("v");
    //model.plotOn(frame, RooFit::Components(poly), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	RooPlot* frame = x.frame();
	RooPlot* framecheck = x.frame();

	if (opt1 == 1) {
		mc.plotOn(frame,RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
	}
	if (opt1 == 2) {
		data.plotOn(frame,RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
	}
	if (opt == 1) {
		pdf.plotOn(frame,RooFit::LineWidth(1));
		pdf.paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
	}

	if (opt == 2) {
		purepdf.plotOn(frame,RooFit::LineWidth(1));
		purepdf.paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		purepdf.plotOn(framecheck,RooFit::Components(exp), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	}

	if (opt == 3) {
		temp_sig.plotOn(frame,RooFit::LineWidth(1));
		temp_sig.paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		temp_sig.plotOn(framecheck,RooFit::Components(histpdf1), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
	} 
	
	RooHist* residuals = frame->residHist("","",true,false); //true = pull, false = center of the bin
	RooPlot* pullFrame = x.frame();
	pullFrame->addPlotable(residuals, "P");
	
	TCanvas *c1 = new TCanvas("c1","",1600,1200);
	c1->Divide(2,2);
	c1->cd(1);
	frame->Draw();
	double meanofBW = bwmean.getVal() - 91.1876;
	double widthofBW = width.getVal() - 2.4955;
	TPaveText *textBox = new TPaveText(0.1, 0.5, 0.3, 0.8, "NDC");
    textBox->SetFillColor(0);
    textBox->SetTextSize(0.03);
    textBox->AddText(Form("Mean value: %.2f", meanofBW));
	textBox->AddText(Form("Width value: %.2f",widthofBW));
    textBox->Draw();
	c1->cd(2);
	pullFrame->Draw();
	c1->cd(3);
	framecheck->Draw();
	

    if (opt == 1 && opt1 == 1 && opt2 == 1) c1->SaveAs("./fitresultplot/Z_mc_cb_bw.pdf");
	if (opt == 1 && opt1 == 1 && opt2 == 2) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_24.pdf");
	if (opt == 1 && opt1 == 1 && opt2 == 3) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_100.pdf");
	if (opt == 1 && opt1 == 2 && opt2 == 1) c1->SaveAs("./fitresultplot/Z_data_cb_bw.pdf");
	if (opt == 1 && opt1 == 2 && opt2 == 2) c1->SaveAs("./fitresultplot/Z_data_cb_bw_24.pdf");
	if (opt == 1 && opt1 == 2 && opt2 == 3) c1->SaveAs("./fitresultplot/Z_data_cb_bw_100.pdf");
	if (opt == 2 && opt1 == 1 && opt2 == 1) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_exp.pdf");
	if (opt == 2 && opt1 == 1 && opt2 == 2) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_exp_24.pdf");
	if (opt == 2 && opt1 == 1 && opt2 == 3) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_exp_100.pdf");
	if (opt == 2 && opt1 == 2 && opt2 == 1) c1->SaveAs("./fitresultplot/Z_data_cb_bw_exp.pdf");
	if (opt == 2 && opt1 == 2 && opt2 == 2) c1->SaveAs("./fitresultplot/Z_data_cb_bw_exp_24.pdf");
	if (opt == 2 && opt1 == 2 && opt2 == 3) c1->SaveAs("./fitresultplot/Z_data_cb_bw_exp_100.pdf");
	if (opt == 3 && opt1 == 1 && opt2 == 1) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_temp.pdf");
	if (opt == 3 && opt1 == 1 && opt2 == 2) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_temp_24.pdf");
	if (opt == 3 && opt1 == 1 && opt2 == 3) c1->SaveAs("./fitresultplot/Z_mc_cb_bw_temp_100.pdf");
	if (opt == 3 && opt1 == 2 && opt2 == 1) c1->SaveAs("./fitresultplot/Z_data_cb_bw_temp.pdf");
	if (opt == 3 && opt1 == 2 && opt2 == 2) c1->SaveAs("./fitresultplot/Z_data_cb_bw_temp_24.pdf");
	if (opt == 3 && opt1 == 2 && opt2 == 3) c1->SaveAs("./fitresultplot/Z_data_cb_bw_temp_100.pdf");
	

}