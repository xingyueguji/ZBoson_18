#include "MC_18.h"
#include "MuonTnP.h"
#include "MCReweight.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"
#include "TEfficiency.h"
#include "ptReweightSpectrum.h"

//C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

void get_ss_bk(){
	TH1::SetDefaultSumw2();
	int nbins = 120;
	TFile *mc_file;
	TFile *data_file;

	mc_file = new TFile("./rootfile/mc_signal.root","UPDATE");
	data_file = new TFile("./rootfile/data_file.root","READ");

	TH1D* mass_array_data_same_sign[11];
	TH1D* mass_array_data_same_sign_with_eff[11];
	TH1D* mass_array_data_same_sign_witheta[11];
	TH1D* mass_array_data_same_sign_witheta_witheff[11];

	TH1D* mass_array_data[11];
	TH1D* mass_array_data_with_eff[11];
	TH1D* mass_array_data_witheta[11];
	TH1D* mass_array_data_witheta_witheff[11];

	TH1D* mass_array[11];
	TH1D* mass_array_with_eff[11];
	TH1D* mass_array_witheta[11];
	TH1D* mass_aray_witheta_witheff[11];

	TH1D* same_sign_ratio[11];
	TH1D* same_sign_ratio_with_eff[11];
	TH1D* same_sign_ratio_witheta[11];
	TH1D* same_sign_ratio_witheta_witheff[11];

	for (int i = 0; i < 11; i++){
		mass_array_data[i] = (TH1D*) data_file->Get(Form("mass_array_data_%i",i));
		mass_array_data_with_eff[i] = (TH1D*) data_file->Get(Form("mass_array_data_with_eff_%i",i));
		mass_array_data_witheta[i] = (TH1D*) data_file->Get(Form("mass_array_data_witheta_%i",i));
		mass_array_data_witheta_witheff[i] = (TH1D*) data_file->Get(Form("mass_array_data_witheta_witheff_%i",i));

		mass_array_data_same_sign[i] = (TH1D*) data_file->Get(Form("mass_array_data_same_sign_%i",i));
		mass_array_data_same_sign_with_eff[i] = (TH1D*) data_file->Get(Form("mass_array_data_same_sign_with_eff_%i",i));
		mass_array_data_same_sign_witheta[i] = (TH1D*) data_file->Get(Form("mass_array_data_same_sign_witheta_%i",i));
		mass_array_data_same_sign_witheta_witheff[i] = (TH1D*) data_file->Get(Form("mass_array_data_same_sign_witheta_witheff_%i",i));

		mass_array[i] = (TH1D*) mc_file->Get(Form("mass_array_%i",i));
		mass_array_with_eff[i] = (TH1D*) mc_file->Get(Form("mass_array_with_eff_%i",i));
		mass_array_witheta[i] = (TH1D*) mc_file->Get(Form("mass_array_witheta_%i",i));
		mass_aray_witheta_witheff[i] = (TH1D*) mc_file->Get(Form("mass_array_witheta_witheff_%i",i));

		same_sign_ratio[i] = new TH1D(Form("same_sign_ratio_%i",i),Form("same_sign_ratio_%i",i),120,60,120);
		same_sign_ratio_with_eff[i] = new TH1D(Form("same_sign_ratio_with_eff_%i",i),Form("same_sign_ratio_with_eff_%i",i),120,60,120);
		same_sign_ratio_witheta[i] = new TH1D(Form("same_sign_ratio_witheta_%i",i),Form("same_sign_ratio_witheta_%i",i),120,60,120);
		same_sign_ratio_witheta_witheff[i] = new TH1D(Form("same_sign_ratio_witheta_witheff_%i",i),Form("same_sign_ratio_witheta_witheff_%i",i),120,60,120);
	}

	for (int i = 0; i < 11; i++){
		for (int j = 1; j<=120; j++){

			double content_mc_os = mass_array[i]->GetBinContent(j);
			double content_data_ss = mass_array_data_same_sign[i]->GetBinContent(j);
			double content_data_os = mass_array_data[i]->GetBinContent(j);
			double ratio;

			double content_mc_os_witheff = mass_array_with_eff[i]->GetBinContent(j);
			double content_data_ss_witheff = mass_array_data_same_sign_with_eff[i]->GetBinContent(j);
			double content_data_os_witheff = mass_array_data_with_eff[i]->GetBinContent(j);
			double ratio_witheff;

			double content_mc_os_witheta = mass_array_witheta[i]->GetBinContent(j);
			double content_data_ss_witheta = mass_array_data_same_sign_witheta[i]->GetBinContent(j);
			double content_data_os_witheta = mass_array_data_witheta[i]->GetBinContent(j);
			double ratio_witheta;

			double content_mc_os_witheta_witheff = mass_aray_witheta_witheff[i]->GetBinContent(j);
			double content_data_ss_witheta_witheff = mass_array_data_same_sign_witheta_witheff[i]->GetBinContent(j);
			double content_data_os_witheta_witheff = mass_array_data_witheta_witheff[i]->GetBinContent(j);
			double ratio_witheta_witheff;

			if (content_data_os >0 ) ratio = content_data_ss / content_data_os;
			else ratio = 0;
			if (content_data_os_witheff >0 ) ratio_witheff = content_data_ss_witheff / content_data_os_witheff;
			else ratio_witheff = 0;
			if (content_data_os_witheta >0 ) ratio_witheta = content_data_ss_witheta / content_data_os_witheta;
			else ratio_witheta = 0;
			if (content_data_os_witheta_witheff >0 ){
				ratio_witheta_witheff = content_data_ss_witheta_witheff / content_data_os_witheta_witheff;
				//cout << "We have content here" << endl;

			} 
			else ratio_witheta_witheff = 0;

			//if (i == 7) cout << "Ratio is " <<  ratio << endl;
			double estimate_mc_ss = content_mc_os * ratio;
			double estimate_mc_ss_witheff = content_mc_os_witheff * ratio_witheff;
			double estimate_mc_ss_witheta = content_mc_os_witheta * ratio_witheta;
			double estimate_mc_ss_witheta_witheff = content_mc_os_witheta_witheff * ratio_witheta_witheff;


			if (estimate_mc_ss >= 0 )same_sign_ratio[i]->SetBinContent(j,estimate_mc_ss);
			if (estimate_mc_ss_witheff >= 0 )same_sign_ratio_with_eff[i]->SetBinContent(j,estimate_mc_ss_witheff);
			if (estimate_mc_ss_witheta >= 0 )same_sign_ratio_witheta[i]->SetBinContent(j,estimate_mc_ss_witheta);
			if (estimate_mc_ss_witheta_witheff >= 0 )same_sign_ratio_witheta_witheff[i]->SetBinContent(j,estimate_mc_ss_witheta_witheff);
			
			//if (i == 7) cout << "Bin content is " <<  estimate_mc_ss << endl;
		}
	}
	mc_file->cd();
	for (int i=0; i<11; i++){
		same_sign_ratio[i]->Write(Form("mc_estimate_ss_%i",i),2);
		same_sign_ratio_with_eff[i]->Write(Form("mc_estimate_ss_with_eff_%i",i),2);
		same_sign_ratio_witheta[i]->Write(Form("mc_estimate_ss_witheta_%i",i),2);
		same_sign_ratio_witheta_witheff[i]->Write(Form("mc_estimate_ss_witheta_witheff_%i",i),2);

	}
	mc_file->Close();
	data_file->Close();

}