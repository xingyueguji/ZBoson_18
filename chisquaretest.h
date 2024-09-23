#include "plotting_helper.h"
class chisquaretest: public plotting_helper{
	public:

	chisquaretest(TString s1, TString s2, TString s3, int type);
	~chisquaretest();
	void calculatechisq(bool isbk);
	void plottingandformatting(int type, bool isbk);
	void bincontentcheck(bool isbk);
	void RebinAll(int x);

	static const int nbins_mass_shift = 21;
	static const int nbins_smear = 21;
	static const int nbins_cent = 11;

	TH1D* h_mc_signal[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D* h_data[nbins_cent];
	TH1D* h_mc_bk[nbins_cent];
	TH1D* h_data_bksub[nbins_cent];
	TH2D* h_chisquare[nbins_cent];

	double chisquarearray[nbins_mass_shift][nbins_smear][nbins_cent] = {};

	TString mcfilepath;
	TString datafilepath;
	TString bkfilepath;

	TCanvas* c_2d_chisquare_ndf[nbins_cent];
	//What I need: invariant mass comparison data vs mc raw, data vs data bk sub.
	TCanvas* c_data_mc_raw[nbins_cent];
	TCanvas* c_data_data_bk[nbins_cent];


	TFile* mcfile;
	TFile* datafile;
	TFile* bkfile;


};

chisquaretest::chisquaretest(TString s1, TString s2, TString s3, int type){

	mcfilepath = s1;
	datafilepath = s2;
	bkfilepath = s3;

	mcfile = new TFile(mcfilepath,"READ");
	datafile = new TFile(datafilepath, "READ");
	bkfile = new TFile(bkfilepath,"READ");

	for (int cent = 0; cent < nbins_cent; cent++){

		c_2d_chisquare_ndf[cent] = new TCanvas(Form("c_2d_chisquare_ndf_%i",cent),"",1600,1200);
		c_data_mc_raw[cent] = new TCanvas(Form("c_data_mc_raw_%i",cent),"",800,600);
		c_data_data_bk[cent] = new TCanvas(Form("c_data_data_bk_%i",cent),"",800,600);

		if (type == 0) h_data[cent] = (TH1D*) datafile ->Get(Form("mass_array_data_with_eff_%i",cent));
		if (type == 1) h_data[cent] = (TH1D*) datafile ->Get(Form("mass_array_data_witheta_witheff_%i",cent));
		h_mc_bk[cent] = (TH1D*) bkfile->Get(Form("Normalized_mc_bk_%i",cent));
		h_chisquare[cent] = new TH2D(Form("h_chisquare_%i",cent),Form("Cent_%i_%i",this->cenlowlimit[cent],this->cenhighlimit[cent]),nbins_mass_shift,-0.5175,0.2175,nbins_smear,-0.0005,0.0205); // Equation here: half bin to the left and half bin to the right, bin width = range / (21-1)

		this->areanormalize(h_data[cent]);

		h_data_bksub[cent] = (TH1D*) h_data[cent]->Clone();
		h_data_bksub[cent]->Add(h_mc_bk[cent],-1);


		for (int shift = 0; shift < nbins_mass_shift; shift++){
			for (int smear = 0; smear < nbins_smear; smear++){
				if (type == 0) h_mc_signal[shift][smear][cent] = (TH1D*) mcfile->Get(Form("modifiedmass_raw_%i_%i_%i",shift,smear,cent));
				if (type == 1) h_mc_signal[shift][smear][cent] = (TH1D*) mcfile->Get(Form("modifiedmass_%i_%i_%i",shift,smear,cent));
				this->areanormalize(h_mc_signal[shift][smear][cent]);
			}
		}
	}
}

void chisquaretest::calculatechisq(bool isbk){
	
	for (int cent = 0; cent < nbins_cent; cent++){
		for (int shift = 0; shift < nbins_mass_shift; shift++){
			for (int smear = 0; smear < nbins_smear; smear ++){
				double chisquarevalue = 0;
				cout << "Now running File " << cent << " mass shift " << shift << " smearing " << smear << endl;
				if (isbk) chisquarevalue = h_data_bksub[cent] -> Chi2Test(h_mc_signal[shift][smear][cent],"WW CHI2/NDF");
				if (!isbk) chisquarevalue = h_data[cent] -> Chi2Test(h_mc_signal[shift][smear][cent],"WW CHI2/NDF");
				//cout << chisquarevalue << endl;
				//chisquarearray[shift][smear][cent] = chisquarevalue;
				std::ostringstream stream;
				stream << std::fixed << std::setprecision(2) << chisquarevalue;
				double formattedBinContent = std::stod(stream.str());
				h_chisquare[cent]->SetBinContent(shift+1,smear+1,formattedBinContent);
			}
		}
	}

}
void chisquaretest::plottingandformatting(int type, bool isbk){

	TString chi2_title;
	TString chi2_saving_path;
	TString data_mc_title;
	TString data_data_title;
	TString data_mc_saving_path;
	TString data_data_saving_path;
	
	if (type == 0){
		if (isbk){
			chi2_title = "Raw_with_bksub_%i_%i";
			data_mc_title = "Raw_with_bksub_%i_%i";
			data_data_title = "Raw_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/raw/bksub/Raw_with_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/raw/bksub/Raw_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_%i_%i.png";
		} 
		if (!isbk){
			chi2_title = "Raw_without_bksub_%i_%i";
			data_mc_title = "Raw_without_bksub_%i_%i";
			data_data_title = "Raw_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/raw/nobksub/Raw_without_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/raw/nobksub/Raw_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_%i_%i.png";
		} 
	}
	if (type == 1){
		if (isbk){
			chi2_title = "Eta_with_bksub_%i_%i";
			data_mc_title = "Eta_with_bksub_%i_%i";
			data_data_title = "Eta_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/eta/bksub/Eta_with_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/eta/bksub/Eta_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_%i_%i.png";
		} 
		if (!isbk){
			chi2_title = "Eta_without_bksub_%i_%i";
			data_mc_title = "Eta_without_bksub_%i_%i";
			data_data_title = "Eta_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/eta/nobksub/Eta_without_bksub_%i_%i.png";
			chi2_saving_path = "./newchi2/chi2plots/eta/nobksub/Eta_with_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/eta/nobksub/Eta_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_%i_%i.png";
		} 
	}

	for (int cent = 0; cent < nbins_cent; cent++){
		//Chi2/ndf plot
		c_2d_chisquare_ndf[cent]->cd();
		gStyle->SetPalette(kRainBow);
		//c_2d_chisquare_ndf[cent]->SetLogz();
		h_chisquare[cent]->SetTitle(Form(chi2_title,this->cenlowlimit[cent],this->cenhighlimit[cent]));
		h_chisquare[cent]->Draw("COLZ");
		h_chisquare[cent]->Draw("TEXTSAME");
		h_chisquare[cent]->GetXaxis()->SetNdivisions(21, 0, 0);
		h_chisquare[cent]->GetYaxis()->SetNdivisions(21, 0, 0);
		h_chisquare[cent]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
    	h_chisquare[cent]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare[cent]->GetXaxis()->SetTitle("Shifted Amount (GeV)");
		h_chisquare[cent]->GetYaxis()->SetTitle("Smeared Amount (Sig)");
		

		for (int j=1; j<=nbins_mass_shift; j++){
			for (int k=1; k<=nbins_smear; k++){
				double xlow = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(k);
        		double xup = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(k);
        		double ylow = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(j);
        		double yup = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(j);

				TBox *box = new TBox(xlow, ylow, xup, yup);
        		box->SetFillStyle(0); // No fill
        		box->SetLineColor(kBlack); // Black border
        		box->SetLineWidth(1); // Border width

        		box->Draw("same");
			}
		}

		Int_t minBinX = -1, minBinY = -1;
		Double_t minContent = h_chisquare[cent]->GetMaximum();
		for (Int_t binX = 1; binX <= h_chisquare[cent]->GetNbinsX(); ++binX) {
        	for (Int_t binY = 1; binY <= h_chisquare[cent]->GetNbinsY(); ++binY) {
            	Double_t content = h_chisquare[cent]->GetBinContent(binX, binY);
            	if (content < minContent) {
              		minContent = content;
              	    minBinX = binX;
                	minBinY = binY;
            	}
        	}
    	}
		
		Double_t xMin = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(minBinX);
   		Double_t xMax = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(minBinX);
    	Double_t yMin = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(minBinY);
    	Double_t yMax = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(minBinY);

		TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
   		box1->SetLineColor(kRed);
    	box1->SetLineWidth(4);
    	box1->SetFillStyle(0);
    	box1->Draw("same");

		c_2d_chisquare_ndf[cent]->SaveAs(Form(chi2_saving_path,this->cenlowlimit[cent],this->cenhighlimit[cent]));

		//This is Data and MC

		c_data_mc_raw[cent]->cd();
		if (isbk) {
			h_data_bksub[cent]->SetTitle(Form(data_mc_title,this->cenlowlimit[cent],this->cenhighlimit[cent]));
			h_data_bksub[cent]->SetMarkerColor(kRed);
			h_data_bksub[cent]->SetMarkerStyle(kFullCircle);
			h_data_bksub[cent]->Draw("P HIST");
			h_mc_signal[14][0][cent]->Draw("P SAME");
		}
		if (!isbk){
			h_data[cent]->SetTitle(Form(data_mc_title,this->cenlowlimit[cent],this->cenhighlimit[cent]));
			h_data[cent]->SetMarkerColor(kRed);
			h_data[cent]->SetMarkerStyle(kFullCircle);
			h_data[cent]->Draw("P HIST");
			h_mc_signal[14][0][cent]->Draw("P SAME");
		}
		c_data_mc_raw[cent]->SaveAs(Form(data_mc_saving_path,this->cenlowlimit[cent],this->cenhighlimit[cent]));


		//This is data vs data_bksub
		c_data_data_bk[cent]->cd();
		c_data_data_bk[cent]->SetLogy();
		h_data[cent]->SetTitle(Form(data_data_title,this->cenlowlimit[cent],this->cenhighlimit[cent]));
		h_data[cent]->SetMarkerColor(kRed);
		h_data[cent]->SetMarkerStyle(kFullCircle);
		h_data[cent]->Draw("P HIST");
		h_data_bksub[cent]->SetMarkerStyle(kFullDotLarge);
		h_data_bksub[cent]->Draw("P SAME");

		c_data_data_bk[cent]->SaveAs(Form(data_data_saving_path,this->cenlowlimit[cent],this->cenhighlimit[cent]));
		
	}


}

void chisquaretest::bincontentcheck(bool isbk){
	int nbins;

	for (int cent = 0; cent < nbins_cent; cent++){
		if (isbk) nbins = h_data[cent]->GetNbinsX();
		if (!isbk) nbins = h_data_bksub[cent]->GetNbinsX();

		for (int bin = 1; bin <= nbins; bin++){
			if (!isbk) {
				if(h_data[cent]->GetBinContent(bin) <= 0) cout << "Warning: " << "Cent " << this->cenlowlimit[cent] << " " << this->cenhighlimit[cent] << " Bin: " << bin << " is <= 0 !" << endl; 
			}
			if (isbk){
				if(h_data_bksub[cent]->GetBinContent(bin) <= 0) cout << "Warning(bksub): " << "Cent " << this->cenlowlimit[cent] << " " << this->cenhighlimit[cent] << " Bin: " << bin << " is <= 0 !" << endl;
			}
		}
	}

}

void chisquaretest::RebinAll(int x){
	for (int cent = 0; cent < nbins_cent; cent++){
		h_data[cent]->Rebin(x);
		h_data_bksub[cent]->Rebin(x);

		for (int shift = 0; shift < nbins_mass_shift; shift++){
			for (int smear = 0; smear < nbins_smear; smear++){
				h_mc_signal[shift][smear][cent]->Rebin(x);
			}
		}

	}

}




