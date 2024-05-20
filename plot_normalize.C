#include "plotting_helper.h"
void plot_normalize(int opt = 1){
	// opt == 1 raw
	// opt == 2 ycut
	// opt == 3 eff
	// opt == 4 y + eff
	// opt == 5 etacut
	// opt == 6 eta + eff

	TFile *mc_signal;
	TFile *data_file;
	TFile *mc_tt;
	TFile *mc_w;

	mc_signal = new TFile("./rootfile/mc_signal.root","UPDATE");
	data_file = new TFile("./rootfile/data_file.root","UPDATE");
	mc_tt = new TFile("./rootfile/mc_tt.root","UPDATE");
	mc_w = new TFile("./rootfile/mc_w.root","UPDATE");

	TFile *rawnormalized = new TFile("./rootfile/normalized/rawfile.root","UPDATE");
	TFile *ycutnormalized = new TFile("./rootfile/normalized/ycutfile.root","UPDATE");
	TFile *effnormalized = new TFile("./rootfile/normalized/efffile.root","UPDATE");
	TFile *ycut_efnormalized = new TFile("./rootfile/normalized/ycut_eff_file.root","UPDATE");
	TFile *etacutnormalized =  new TFile("./rootfile/normalized/etacut_file.root","UPDATE");
	TFile *etacut_effnormalized = new TFile("./rootfile/normalized/etacut_eff_file.root","UPDATE");


	TH1D *mass_array_W[11];
	TH1D *mass_array_tt[11];
	TH1D *mass_array_tau[11];
	TH1D *mass_array_signal[11];
	TH1D *mass_array_data[11];
	TH1D *mass_array_samesign[11];

	for (int i=0; i<11; i++){
		if (opt == 1) {

			mass_array_W[i] = (TH1D*) mc_w->Get(Form("mass_array_%i",i));
			mass_array_tt[i] = (TH1D*) mc_tt->Get(Form("mass_array_%i",i));
			mass_array_tau[i] = (TH1D*) mc_signal->Get(Form("mass_array_tau_%i",i));
			mass_array_signal[i] = (TH1D*) mc_signal->Get(Form("mass_array_%i",i));
			mass_array_data[i] = (TH1D*) data_file->Get(Form("mass_array_data_%i",i));
			mass_array_samesign[i] = (TH1D*) mc_signal->Get(Form("mc_estimate_ss_%i",i));

		}
		if (opt == 2) {

			mass_array_W[i] = (TH1D*) mc_w->Get(Form("mass_array_withy_%i",i));
			mass_array_tt[i] = (TH1D*) mc_tt->Get(Form("mass_array_withy_%i",i));
			mass_array_tau[i] = (TH1D*) mc_signal->Get(Form("mass_array_tau_withy_%i",i));
			mass_array_signal[i] = (TH1D*) mc_signal->Get(Form("mass_array_withy_%i",i));
			mass_array_data[i] = (TH1D*) data_file->Get(Form("mass_array_data_withy_%i",i));
			mass_array_samesign[i] = (TH1D*) mc_signal->Get(Form("mc_estimate_ss_withy_%i",i));

		}
		if (opt == 3) {

			mass_array_W[i] = (TH1D*) mc_w->Get(Form("mass_array_with_eff_%i",i));
			mass_array_tt[i] = (TH1D*) mc_tt->Get(Form("mass_array_with_eff_%i",i));
			mass_array_tau[i] = (TH1D*) mc_signal->Get(Form("mass_array_tau_with_eff_%i",i));
			mass_array_signal[i] = (TH1D*) mc_signal->Get(Form("mass_array_with_eff_%i",i));
			mass_array_data[i] = (TH1D*) data_file->Get(Form("mass_array_data_with_eff_%i",i));
			mass_array_samesign[i] = (TH1D*) mc_signal->Get(Form("mc_estimate_ss_with_eff_%i",i));

		}
		if (opt == 4) {

			mass_array_W[i] = (TH1D*) mc_w->Get(Form("mass_array_withy_witheff_%i",i));
			mass_array_tt[i] = (TH1D*) mc_tt->Get(Form("mass_array_withy_witheff_%i",i));
			mass_array_tau[i] = (TH1D*) mc_signal->Get(Form("mass_array_tau_withy_witheff_%i",i));
			mass_array_signal[i] = (TH1D*) mc_signal->Get(Form("mass_array_withy_witheff_%i",i));
			mass_array_data[i] = (TH1D*) data_file->Get(Form("mass_array_data_withy_witheff_%i",i));
			mass_array_samesign[i] = (TH1D*) mc_signal->Get(Form("mc_estimate_ss_withy_witheff_%i",i));

		}
		if (opt == 5){

			mass_array_W[i] = (TH1D*) mc_w->Get(Form("mass_array_witheta_%i",i));
			mass_array_tt[i] = (TH1D*) mc_tt->Get(Form("mass_array_witheta_%i",i));
			mass_array_tau[i] = (TH1D*) mc_signal->Get(Form("mass_array_tau_witheta_%i",i));
			mass_array_signal[i] = (TH1D*) mc_signal->Get(Form("mass_array_witheta_%i",i));
			mass_array_data[i] = (TH1D*) data_file->Get(Form("mass_array_data_witheta_%i",i));
			mass_array_samesign[i] = (TH1D*) mc_signal->Get(Form("mc_estimate_ss_witheta_%i",i));

		}
		if (opt == 6){

			mass_array_W[i] = (TH1D*) mc_w->Get(Form("mass_array_witheta_witheff_%i",i));
			mass_array_tt[i] = (TH1D*) mc_tt->Get(Form("mass_array_witheta_witheff_%i",i));
			mass_array_tau[i] = (TH1D*) mc_signal->Get(Form("mass_array_tau_witheta_witheff_%i",i));
			mass_array_signal[i] = (TH1D*) mc_signal->Get(Form("mass_array_witheta_witheff_%i",i));
			mass_array_data[i] = (TH1D*) data_file->Get(Form("mass_array_data_witheta_witheff_%i",i));
			mass_array_samesign[i] = (TH1D*) mc_signal->Get(Form("mc_estimate_ss_witheta_witheff_%i",i));

		}
	}

	TH1D* h_DY_weight;
	TH1D* h_tt_weight;
	TH1D* h_W_weight;

	h_DY_weight = (TH1D*) mc_signal->Get("event_weight");
	h_tt_weight = (TH1D*) mc_tt->Get("event_weight");
	h_W_weight = (TH1D*) mc_w->Get("event_weight");

	double w_norm = h_W_weight->Integral();
	double tt_norm = h_tt_weight->Integral();
	double dy_norm = h_DY_weight->Integral();

	plotting_helper *ovo = new plotting_helper();
	ovo->setTDRStyle();

	TH1D* h_normalized_data[11];
	TH1D* h_normalized_mc[11];
	TH1D* h_normalized_mc_bk[11];
	TH1D* h_total_mc_30bins[11];

	for (int i=0; i<11; i++){
		ovo->luminormalize(mass_array_W[i],2,w_norm);
		ovo->luminormalize(mass_array_tt[i],3,tt_norm);
		ovo->luminormalize(mass_array_tau[i],1,dy_norm);
		ovo->luminormalize(mass_array_signal[i],1,dy_norm);
		ovo->luminormalize(mass_array_samesign[i],1,dy_norm);

		h_normalized_data[i] = (TH1D*) mass_array_data[i]->Clone();

		h_normalized_mc[i] = (TH1D*) mass_array_signal[i]->Clone();
		h_normalized_mc[i]->Add(mass_array_W[i],1);
		h_normalized_mc[i]->Add(mass_array_tt[i],1);
		h_normalized_mc[i]->Add(mass_array_tau[i],1);
		h_normalized_mc[i]->Add(mass_array_samesign[i],1);

		h_normalized_mc_bk[i] = (TH1D*) mass_array_W[i] ->Clone();
		h_normalized_mc_bk[i]->Add(mass_array_tt[i],1);
		h_normalized_mc_bk[i]->Add(mass_array_tau[i],1);
		h_normalized_mc_bk[i]->Add(mass_array_samesign[i],1);

		ovo->areanormalize(h_normalized_data[i]);
		ovo->areanormalize(h_normalized_mc[i]);
		ovo->areanormalize(h_normalized_mc_bk[i]);

		mass_array_W[i]->Rebin(4);
		mass_array_tt[i]->Rebin(4);
		mass_array_tau[i]->Rebin(4);
		mass_array_signal[i]->Rebin(4);
		mass_array_data[i]->Rebin(4);
		mass_array_samesign[i]->Rebin(4);

		h_total_mc_30bins[i] = (TH1D*) mass_array_W[i]->Clone();
		h_total_mc_30bins[i]->Add(mass_array_tt[i]);
		h_total_mc_30bins[i]->Add(mass_array_tau[i]);
		h_total_mc_30bins[i]->Add(mass_array_signal[i]);
		h_total_mc_30bins[i]->Add(mass_array_samesign[i]);

		double normalization_factor_mc_30bins = h_total_mc_30bins[i]->Integral("width");
		mass_array_W[i]->Scale(1/normalization_factor_mc_30bins);
		mass_array_tt[i]->Scale(1/normalization_factor_mc_30bins);
		mass_array_tau[i]->Scale(1/normalization_factor_mc_30bins);
		mass_array_signal[i]->Scale(1/normalization_factor_mc_30bins);
		mass_array_samesign[i]->Scale(1/normalization_factor_mc_30bins);

		ovo->areanormalize(mass_array_data[i]);

		ovo->compositeplot(mass_array_data[i],mass_array_signal[i],mass_array_samesign[i],mass_array_tau[i],mass_array_W[i],mass_array_tt[i],i,opt);
		if (opt == 1)ovo->savehistogram(h_normalized_mc[i],h_normalized_data[i],h_normalized_mc_bk[i],i,rawnormalized);
		if (opt == 2)ovo->savehistogram(h_normalized_mc[i],h_normalized_data[i],h_normalized_mc_bk[i],i,ycutnormalized);
		if (opt == 3)ovo->savehistogram(h_normalized_mc[i],h_normalized_data[i],h_normalized_mc_bk[i],i,effnormalized);
		if (opt == 4)ovo->savehistogram(h_normalized_mc[i],h_normalized_data[i],h_normalized_mc_bk[i],i,ycut_efnormalized);
		if (opt == 5)ovo->savehistogram(h_normalized_mc[i],h_normalized_data[i],h_normalized_mc_bk[i],i,etacutnormalized);
		if (opt == 6)ovo->savehistogram(h_normalized_mc[i],h_normalized_data[i],h_normalized_mc_bk[i],i,etacut_effnormalized);
		
	}


}