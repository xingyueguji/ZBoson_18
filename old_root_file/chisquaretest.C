#include "plotting_helper.h"
void chisquaretest(int type = 1,bool isbksub = 1){
	//type == 0 raw
	//type == 1 eta cut

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	//Let's shift the mc

	const int nbins_mass_shift = 21;
	const int nbins_smear = 21;
	const int nbins_cent = 11;

	TH1D* mc_hist[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D* data_hist[11];
	TH1D* data_hist_bksub[11];
	TH1D* bk_hist[11];
	double chisquarearray[nbins_mass_shift][nbins_smear][nbins_cent] = {};

	TH2D* h_chisquare[nbins_cent];
	TCanvas *c_chisquare[nbins_cent];
	TCanvas *c_bksubcheck[nbins_cent];

	string Filenames[nbins_cent];

	plotting_helper *ovo = new plotting_helper();
	//ovo->setTDRStyle();

	TFile *f1 = new TFile("./rootfile/modified_signal.root","READ");
	TFile *f2 = new TFile("./rootfile/data_file.root","READ");
	TFile *f3;
	if (type == 0) f3 = new TFile("./rootfile/normalized/efffile.root");
	if (type == 1) f3 = new TFile("./rootfile/normalized/etacut_eff_file.root");

	for (int i = 0; i < 11; i++) {
		if (type == 0) {
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_with_eff_%i",i));
			bk_hist[i] = (TH1D*) f3->Get(Form("Normalized_mc_bk_%i",i));
		}
		if (type == 1) {
			data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_witheta_witheff_%i",i));
			bk_hist[i] = (TH1D*) f3->Get(Form("Normalized_mc_bk_%i",i));
		}
		
		ovo->areanormalize(data_hist[i]);
		data_hist_bksub[i] = (TH1D*) data_hist[i]->Clone();
		
		h_chisquare[i] = new TH2D(Form("h_chisquare_%i",i),Form("Cent_%i_%i",ovo->cenlowlimit[i],ovo->cenhighlimit[i]),nbins_mass_shift,-0.5175,0.2175,nbins_smear,-0.0005,0.0205); // Equation here: half bin to the left and half bin to the right, bin width = range / (21-1)
		c_chisquare[i] = new TCanvas(Form("c_%i",i),"",1600,1200);
		c_bksubcheck[i] = new TCanvas(Form("c_1_%i",i),"",1600,1200);
		c_chisquare[i]->SetLogz();

		c_bksubcheck[i]->SetLogy();
		c_bksubcheck[i]->cd();

		data_hist[i]->SetMarkerStyle(20);
		data_hist[i]->Draw("P");
		data_hist_bksub[i]->Add(bk_hist[i],-1);
		data_hist_bksub[i]->SetLineColor(kRed);
		data_hist_bksub[i]->SetMarkerStyle(24);
		data_hist_bksub[i]->Draw("P SAME");

		if (isbksub) data_hist[i]->Add(bk_hist[i],-1);

	}

	for (int i = 0; i < nbins_mass_shift; i++){
		for (int j = 0; j < nbins_smear; j++){
			for (int k = 0; k < nbins_cent; k++){
				if (type == 0) mc_hist[i][j][k] = (TH1D*) f1->Get(Form("modifiedmass_raw_%i_%i_%i",i,j,k));
				if (type == 1) mc_hist[i][j][k] = (TH1D*) f1->Get(Form("modifiedmass_%i_%i_%i",i,j,k));
				ovo->areanormalize(mc_hist[i][j][k]);
			}
		}
	}
	

	for (int i = 0; i < nbins_mass_shift; i++){
		for (int j = 0; j < nbins_smear; j++){
			for (int k = 0; k < nbins_cent; k++){
				cout << "Now running File " << k << " mass shift " << i << " smearing " << j << endl;
				double chi2value = 0;
				chi2value = data_hist[k] -> Chi2Test(mc_hist[i][j][k],"WW CHI2/NDF");

				cout << chi2value << endl;
				chisquarearray[i][j][k] = chi2value;
			}
		}
	}

	for (int i=0; i< nbins_cent; i++){
		for (int j=0; j<nbins_smear; j++){
			for (int k=0; k<nbins_mass_shift; k++){

				std::ostringstream stream;
				stream << std::fixed << std::setprecision(2) << chisquarearray[k][j][i]; // Change precision as needed
				double formattedBinContent = std::stod(stream.str());

				h_chisquare[i]->SetBinContent(k+1,j+1,formattedBinContent);
				//outFiles[i] << chisquarearray[k][j][i];
				//if (k < 20) outFiles[i] << " ";
			}
			//if (j < 20) outFiles[i] << endl;
		}
	}

	for (int i=0; i < nbins_cent; i++){
		c_chisquare[i]->cd();
		if (type == 0) {
			if (isbksub) h_chisquare[i]->SetTitle(Form("raw_cent_%i_%i_bksub",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
			if (!isbksub) h_chisquare[i]->SetTitle(Form("raw_cent_%i_%i",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
		}
		if (type == 1) {
			if (isbksub) h_chisquare[i]->SetTitle(Form("etacut_cent_%i_%i_bksub",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
			if (!isbksub) h_chisquare[i]->SetTitle(Form("etacut_cent_%i_%i",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
		}
		h_chisquare[i]->Draw("COLZ");
		h_chisquare[i]->Draw("TEXTSAME");

		h_chisquare[i]->GetXaxis()->SetNdivisions(21, 0, 0);
		h_chisquare[i]->GetYaxis()->SetNdivisions(21, 0, 0);
		h_chisquare[i]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
    	h_chisquare[i]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare[i]->GetXaxis()->SetTitle("Shifted Amount (GeV)");
		h_chisquare[i]->GetYaxis()->SetTitle("Smeared Amount (Sig)");



		for (int j=1; j<=nbins_mass_shift; j++){
			for (int k=1; k<=nbins_smear; k++){
				double xlow = h_chisquare[i]->GetXaxis()->GetBinLowEdge(k);
        		double xup = h_chisquare[i]->GetXaxis()->GetBinUpEdge(k);
        		double ylow = h_chisquare[i]->GetYaxis()->GetBinLowEdge(j);
        		double yup = h_chisquare[i]->GetYaxis()->GetBinUpEdge(j);

				TBox *box = new TBox(xlow, ylow, xup, yup);
        		box->SetFillStyle(0); // No fill
        		box->SetLineColor(kBlack); // Black border
        		box->SetLineWidth(1); // Border width

        		box->Draw("same");
			}
		}
		
		Int_t minBinX = -1, minBinY = -1;
		Double_t minContent = h_chisquare[i]->GetMaximum();
		for (Int_t binX = 1; binX <= h_chisquare[i]->GetNbinsX(); ++binX) {
        	for (Int_t binY = 1; binY <= h_chisquare[i]->GetNbinsY(); ++binY) {
            	Double_t content = h_chisquare[i]->GetBinContent(binX, binY);
            	if (content < minContent) {
              		minContent = content;
              	    minBinX = binX;
                	minBinY = binY;
            	}
        	}
    	}
		
		Double_t xMin = h_chisquare[i]->GetXaxis()->GetBinLowEdge(minBinX);
   		Double_t xMax = h_chisquare[i]->GetXaxis()->GetBinUpEdge(minBinX);
    	Double_t yMin = h_chisquare[i]->GetYaxis()->GetBinLowEdge(minBinY);
    	Double_t yMax = h_chisquare[i]->GetYaxis()->GetBinUpEdge(minBinY);

		TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
   		box1->SetLineColor(kRed);
    	box1->SetLineWidth(4);
    	box1->SetFillStyle(0);
    	box1->Draw("same");

		if (type == 0 && isbksub){
			c_chisquare[i]->SaveAs(Form("./chisquaretest/plots/bksub/raw/Cent_%i_%i_bksub.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
			c_bksubcheck[i]->SaveAs(Form("./chisquaretest/sancheckplots/bksub/raw/Cent_%i_%i_bksub.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
		}
		if (type == 0 && !isbksub){
			c_chisquare[i]->SaveAs(Form("./chisquaretest/plots/raw/Cent_%i_%i.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
			c_bksubcheck[i]->SaveAs(Form("./chisquaretest/sancheckplots/raw/Cent_%i_%i.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
		} 
		if (type == 1 && isbksub){
			c_chisquare[i]->SaveAs(Form("./chisquaretest/plots/bksub/etacut/Cent_%i_%i_bksub.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
			c_bksubcheck[i]->SaveAs(Form("./chisquaretest/sancheckplots/bksub/eta/Cent_%i_%i_bksub.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
		}
		if (type == 1 && !isbksub){
			c_chisquare[i]->SaveAs(Form("./chisquaretest/plots/etacut/Cent_%i_%i.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
			c_bksubcheck[i]->SaveAs(Form("./chisquaretest/sancheckplots/eta/Cent_%i_%i.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
		}
		
	}
}