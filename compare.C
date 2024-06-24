#include "plotting_helper.h"
void compare(int type = 1){
	// type 1 == raw
	// type 2 == raw + eff
	// type 3 == y
	// type 4 == y + eff
	// type 5 == eta
	// type 6 == eta + eff
	TString typestring;
	if (type == 1) typestring = "raw";
	if (type == 2) typestring = "raw+eff";
	if (type == 3) typestring = "y";
	if (type == 4) typestring = "y+eff";
	if (type == 5) typestring = "eta";
	if (type == 6) typestring = "eta+eff";


	plotting_helper *ovo = new plotting_helper();
	ovo->setTDRStyle();

	// I am going to add data + mc histogram and data + mc roodataset

	// Here's histogram
	TH1D* mc_hist[11];
	TH1D* data_hist[11];
	TH1D* diff_hist[11];
	TH1* normalized_data_dataset[11];
	TH1* normalized_mc_dataset[11];
	RooDataSet* data_dataset[11];
	RooDataSet* mc_dataset[11];

	TFile *f1 = new TFile("./rootfile/mc_signal.root","READ");
	TFile *f2 = new TFile("./rootfile/data_file.root","READ");

	if (f1->IsOpen() && !f1->IsZombie()) {
        std::cout << "F1 File successfully opened." << std::endl;
    } else {
        std::cerr << "Error opening file F1." << std::endl;
        // Handle the error (e.g., cleanup, exit, etc.)
    }
	if (f2->IsOpen() && !f2->IsZombie()) {
        std::cout << "F2 File successfully opened." << std::endl;
    } else {
        std::cerr << "Error opening file F2." << std::endl;
        // Handle the error (e.g., cleanup, exit, etc.)
    }



	TCanvas *c1[11];
	TCanvas *c2[11];
	TCanvas *c3[11];
	TCanvas *c4[11];

	RooRealVar *x = new RooRealVar("roomass","roomass",60,120);
	RooPlot* frame[11];
	RooPlot* frame_1[11];
	
	TString savingpathdata = "./manualplot/data/";
	TString savingpathmc = "./manualplot/mc/";
	TString savingpathoverlay = "./manualplot/overlay/";
	TString savingpathdiff = "./manualplot/diff/";
	TString filename;

	TLegend* legend[11];

	for (int i = 0; i < 11; i++) {

		c1[i] = new TCanvas(Form("c1_%i",i),"",800,600);
		c2[i] = new TCanvas(Form("c2_%i",i),"",800,600);
		c3[i] = new TCanvas(Form("c3_%i",i),"",800,600);
		c4[i] = new TCanvas(Form("c4_%i",i),"",800,600);
		frame[i] = x->frame();
		frame_1[i] = x->frame();
		legend[i] = new TLegend(0.75, 0.75, 0.9, 0.9);

		if (type == 1){
			mc_hist[i] = (TH1D*) f1->Get(Form("mass_array_%i",i));
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_%i",i));
			data_dataset[i] = (RooDataSet*) f2->Get(Form("roodata_%i",i));
			mc_dataset[i] = (RooDataSet*) f1->Get(Form("rooreco_%i",i));
			filename = "raw_%i_%i.png";
		}
		if (type == 2){
			mc_hist[i] = (TH1D*) f1->Get(Form("mass_array_with_eff_%i",i));
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_with_eff_%i",i));
			data_dataset[i] = (RooDataSet*) f2->Get(Form("roodata_eff_%i",i));
			mc_dataset[i] = (RooDataSet*) f1->Get(Form("rooreco_eff_%i",i));
			filename = "raw+eff_%i_%i.png";
		}
		if (type == 3){
			mc_hist[i] = (TH1D*) f1->Get(Form("mass_array_withy_%i",i));
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_withy_%i",i));
			data_dataset[i] = (RooDataSet*) f2->Get(Form("roodata_y_%i",i));
			mc_dataset[i] = (RooDataSet*) f1->Get(Form("rooreco_y_%i",i));
			filename = "y_%i_%i.png";
		}
		if (type == 4){
			mc_hist[i] = (TH1D*) f1->Get(Form("mass_array_withy_witheff_%i",i));
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_withy_witheff_%i",i));
			data_dataset[i] = (RooDataSet*) f2->Get(Form("roodata_y_eff_%i",i));
			mc_dataset[i] = (RooDataSet*) f1->Get(Form("rooreco_y_eff_%i",i));
			filename = "y+eff_%i_%i.png";
		}
		if (type == 5){
			mc_hist[i] = (TH1D*) f1->Get(Form("mass_array_witheta_%i",i));
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_witheta_%i",i));
			data_dataset[i] = (RooDataSet*) f2->Get(Form("roodata_eta_%i",i));
			mc_dataset[i] = (RooDataSet*) f1->Get(Form("rooreco_eta_%i",i));
			filename = "eta_%i_%i.png";
		}
		if (type == 6){
			mc_hist[i] = (TH1D*) f1->Get(Form("mass_array_witheta_witheff_%i",i));
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_witheta_witheff_%i",i));
			data_dataset[i] = (RooDataSet*) f2->Get(Form("roodata_eta_eff_%i",i));
			mc_dataset[i] = (RooDataSet*) f1->Get(Form("rooreco_eta_eff_%i",i));
			filename = "eta+eff_%i_%i.png";
		}
	}

	// First check if roodataset and histograms are the same.
	// Normalize all first.
	// RooDataset use plotOn(frame, RooFit::Normalization(1.0, RooAbsReal::Relative)); so it normalized to unit area.

	for (int i=0; i<11; i++){
		ovo->areanormalize(mc_hist[i]);
		ovo->areanormalize(data_hist[i]);

		normalized_data_dataset[i] = (TH1*) data_dataset[i]->createHistogram(Form("data_%i",i),*x, RooFit::Binning(120));
		normalized_mc_dataset[i] = (TH1*) mc_dataset[i]->createHistogram(Form("mc_%i",i),*x, RooFit::Binning(120));

		ovo->areanormalize_TH1(normalized_data_dataset[i]);
		ovo->areanormalize_TH1(normalized_mc_dataset[i]);
	}

	for (int i=0; i<11; i++){
		c1[i]->cd();

		mc_hist[i]->SetLineColor(kRed);
		normalized_mc_dataset[i]->SetMarkerStyle(21);
		normalized_mc_dataset[i]->Draw("P");
		mc_hist[i]->Draw("SAME");

		c2[i]->cd();

		data_hist[i]->SetLineColor(kRed);
		normalized_data_dataset[i]->SetMarkerStyle(21);
		normalized_data_dataset[i]->Draw("P");
		data_hist[i]->Draw("SAME");

		c3[i]->cd();

		data_hist[i]->SetLineColor(kBlack);
		data_hist[i]->Draw();
		mc_hist[i]->Draw("SAME");

		legend[i]->AddEntry(data_hist[i],"DATA " + typestring, "lep");
		legend[i]->AddEntry(mc_hist[i],"MC " + typestring, "lep");

		legend[i]->Draw();

		c4[i]->cd();
		diff_hist[i] = (TH1D*) data_hist[i]->Clone();
		diff_hist[i]->Add(mc_hist[i],-1);
		diff_hist[i]->Draw();

	}

	for (int i=0; i<11; i++){
		TString fullpathdata = savingpathdata + Form(filename,ovo->cenlowlimit[i],ovo->cenhighlimit[i]);
		TString fullpathmc = savingpathmc + Form(filename,ovo->cenlowlimit[i],ovo->cenhighlimit[i]);
		TString fullpathoverlay = savingpathoverlay + Form(filename,ovo->cenlowlimit[i],ovo->cenhighlimit[i]);
		TString fullpathdiff = savingpathdiff + Form(filename,ovo->cenlowlimit[i],ovo->cenhighlimit[i]);
		c1[i]->SaveAs(fullpathmc);
		c2[i]->SaveAs(fullpathdata);
		c3[i]->SaveAs(fullpathoverlay);
		c4[i]->SaveAs(fullpathdiff);
	}




	








}