#include "plotting_helper.h"
void chisquaretest(){

	TH1::SetDefaultSumw2();

	//Let's shift the mc

	const int nbins_mass_shift = 21;
	const int nbins_smear = 21;
	const int nbins_cent = 11;

	TH1D* mc_hist[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D* data_hist[11];
	double chisquarearray[nbins_mass_shift][nbins_smear][nbins_cent] = {};

	TH2D* h_chisquare[nbins_cent];
	TCanvas *c_chisquare[nbins_cent];

	ofstream outFiles[nbins_cent];
	string Filenames[nbins_cent];

	plotting_helper *ovo = new plotting_helper();
	//ovo->setTDRStyle();

	TFile *f1 = new TFile("./rootfile/modified_signal.root","READ");
	TFile *f2 = new TFile("./rootfile/data_file.root","READ");

	for (int i = 0; i < 11; i++) {
		data_hist[i] = (TH1D*) f2->Get(Form("mass_array_data_witheta_witheff_%i",i));
		Filenames[i] = Form("./chisquaretest/Cent_%i_%i.txt",ovo->cenlowlimit[i],ovo->cenhighlimit[i]);
		outFiles[i].open(Filenames[i]);
		h_chisquare[i] = new TH2D(Form("h_chisquare_%i",i),Form("Cent_%i_%i",ovo->cenlowlimit[i],ovo->cenhighlimit[i]),nbins_mass_shift,-0.45,0.55,nbins_smear,-0.05,1.05);
		c_chisquare[i] = new TCanvas(Form("c_%i",i),"",800,600);
	}

	for (int i = 0; i < nbins_mass_shift; i++){
		for (int j = 0; j < nbins_smear; j++){
			for (int k = 0; k < nbins_cent; k++){
				mc_hist[i][j][k] = (TH1D*) f1->Get(Form("modifiedmass_%i_%i_%i",i,j,k));
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
				h_chisquare[i]->SetBinContent(k+1,j+1,chisquarearray[k][j][i]);
				outFiles[i] << chisquarearray[k][j][i];
				if (k < 20) outFiles[i] << " ";
			}
			if (j < 20) outFiles[i] << endl;
		}
	}

	for (int i=0; i < nbins_cent; i++){
		c_chisquare[i]->cd();
		h_chisquare[i]->Draw("COLZ");
		c_chisquare[i]->SaveAs(Form("./chisquaretest/plots/Cent_%i_%i.png",ovo->cenlowlimit[i],ovo->cenhighlimit[i]));
		outFiles[i].close();
	}
}