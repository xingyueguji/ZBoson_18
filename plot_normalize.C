#include "plotting_helper.h"
void plot_normalize(int opt = 7)
{
	// opt == 1 nominal
	// opt == 2 tnpU
	// opt == 3 tnpD
	// opt == 4 Acooff
	// opt == 6 Mass range

	// now raw and eta are saved into same file, only differ by version

	TFile *mc_signal;
	TFile *data_file;
	TFile *mc_tt;
	TFile *mc_w;

	mc_signal = new TFile("./rootfile/mc_signal.root", "READ");
	data_file = new TFile("./rootfile/data_file.root", "READ");
	mc_tt = new TFile("./rootfile/mc_tt.root", "READ");
	mc_w = new TFile("./rootfile/mc_w.root", "READ");

	TFile *f_FA_nominal = new TFile("./rootfile/normalized/FA_nominal.root", "UPDATE");
	TFile *f_FA_tnpU = new TFile("./rootfile/normalized/FA_tnpU.root", "UPDATE");
	TFile *f_FA_tnpD = new TFile("./rootfile/normalized/FA_tnpD.root", "UPDATE");
	TFile *f_FA_Acooff = new TFile("./rootfile/normalized/FA_acooff.root", "UPDATE");
	TFile *f_FA_mass_range = new TFile("./rootfile/normalized/FA_mass_range.root", "UPDATE");

	TFile *f_Eta_nominal = new TFile("./rootfile/normalized/Eta_nominal.root", "UPDATE");
	TFile *f_Eta_tnpU = new TFile("./rootfile/normalized/Eta_tnpU.root", "UPDATE");
	TFile *f_Eta_tnpD = new TFile("./rootfile/normalized/Eta_tnpD.root", "UPDATE");
	TFile *f_Eta_Acooff = new TFile("./rootfile/normalized/Eta_acooff.root", "UPDATE");
	TFile *f_Eta_mass_range = new TFile("./rootfile/normalized/Eta_mass_range.root", "UPDATE");

	TH1D *h_FA_mcsignal[11];
	TH1D *h_FA_mcW[11];
	TH1D *h_FA_mctt[11];
	TH1D *h_FA_mctau[11];
	TH1D *h_FA_data[11];
	TH1D *h_FA_samesign[11];

	TH1D *h_Eta_mcsignal[11];
	TH1D *h_Eta_mcW[11];
	TH1D *h_Eta_mctt[11];
	TH1D *h_Eta_mctau[11];
	TH1D *h_Eta_data[11];
	TH1D *h_Eta_samesign[11];

	for (int i = 0; i < 11; i++)
	{
		if (opt == 1)
		{
			h_FA_mcsignal[i] = (TH1D *)mc_signal->Get(Form("FA_nominal_%i", i));
			h_FA_mcW[i] = (TH1D *)mc_w->Get(Form("FA_nominal_%i", i));
			h_FA_mctt[i] = (TH1D *)mc_tt->Get(Form("FA_nominal_%i", i));
			h_FA_mctau[i] = (TH1D *)mc_signal->Get(Form("FA_tau_nominal_%i", i));
			h_FA_data[i] = (TH1D *)data_file->Get(Form("FA_nominal_%i", i));
			h_FA_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_FA_nominal_%i", i));

			h_Eta_mcsignal[i] = (TH1D *)mc_signal->Get(Form("Eta_nominal_%i", i));
			h_Eta_mcW[i] = (TH1D *)mc_w->Get(Form("Eta_nominal_%i", i));
			h_Eta_mctt[i] = (TH1D *)mc_tt->Get(Form("Eta_nominal_%i", i));
			h_Eta_mctau[i] = (TH1D *)mc_signal->Get(Form("Eta_tau_nominal_%i", i));
			h_Eta_data[i] = (TH1D *)data_file->Get(Form("Eta_nominal_%i", i));
			h_Eta_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_Eta_nominal_%i", i));
		}
		if (opt == 2)
		{

			h_FA_mcsignal[i] = (TH1D *)mc_signal->Get(Form("FA_nominal_%i", i));
			h_FA_mcW[i] = (TH1D *)mc_w->Get(Form("FA_nominal_%i", i));
			h_FA_mctt[i] = (TH1D *)mc_tt->Get(Form("FA_nominal_%i", i));
			h_FA_mctau[i] = (TH1D *)mc_signal->Get(Form("FA_tau_nominal_%i", i));
			h_FA_data[i] = (TH1D *)data_file->Get(Form("FA_tnpU_%i", i));
			h_FA_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_FA_tnpU_%i", i));

			h_Eta_mcsignal[i] = (TH1D *)mc_signal->Get(Form("Eta_nominal_%i", i));
			h_Eta_mcW[i] = (TH1D *)mc_w->Get(Form("Eta_nominal_%i", i));
			h_Eta_mctt[i] = (TH1D *)mc_tt->Get(Form("Eta_nominal_%i", i));
			h_Eta_mctau[i] = (TH1D *)mc_signal->Get(Form("Eta_tau_nominal_%i", i));
			h_Eta_data[i] = (TH1D *)data_file->Get(Form("Eta_tnpU_%i", i));
			h_Eta_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_Eta_tnpU_%i", i));
		}
		if (opt == 3)
		{

			h_FA_mcsignal[i] = (TH1D *)mc_signal->Get(Form("FA_nominal_%i", i));
			h_FA_mcW[i] = (TH1D *)mc_w->Get(Form("FA_nominal_%i", i));
			h_FA_mctt[i] = (TH1D *)mc_tt->Get(Form("FA_nominal_%i", i));
			h_FA_mctau[i] = (TH1D *)mc_signal->Get(Form("FA_tau_nominal_%i", i));
			h_FA_data[i] = (TH1D *)data_file->Get(Form("FA_tnpD_%i", i));
			h_FA_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_FA_tnpD_%i", i));

			h_Eta_mcsignal[i] = (TH1D *)mc_signal->Get(Form("Eta_nominal_%i", i));
			h_Eta_mcW[i] = (TH1D *)mc_w->Get(Form("Eta_nominal_%i", i));
			h_Eta_mctt[i] = (TH1D *)mc_tt->Get(Form("Eta_nominal_%i", i));
			h_Eta_mctau[i] = (TH1D *)mc_signal->Get(Form("Eta_tau_nominal_%i", i));
			h_Eta_data[i] = (TH1D *)data_file->Get(Form("Eta_tnpD_%i", i));
			h_Eta_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_Eta_tnpD_%i", i));
		}
		if (opt == 4)
		{

			h_FA_mcsignal[i] = (TH1D *)mc_signal->Get(Form("FA_AcoOff_%i", i));
			h_FA_mcW[i] = (TH1D *)mc_w->Get(Form("FA_AcoOff_%i", i));
			h_FA_mctt[i] = (TH1D *)mc_tt->Get(Form("FA_AcoOff_%i", i));
			h_FA_mctau[i] = (TH1D *)mc_signal->Get(Form("FA_tau_AcoOff_%i", i));
			h_FA_data[i] = (TH1D *)data_file->Get(Form("FA_AcoOff_%i", i));
			h_FA_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_FA_acooff_%i", i));

			h_Eta_mcsignal[i] = (TH1D *)mc_signal->Get(Form("Eta_AcoOff_%i", i));
			h_Eta_mcW[i] = (TH1D *)mc_w->Get(Form("Eta_AcoOff_%i", i));
			h_Eta_mctt[i] = (TH1D *)mc_tt->Get(Form("Eta_AcoOff_%i", i));
			h_Eta_mctau[i] = (TH1D *)mc_signal->Get(Form("Eta_tau_AcoOff_%i", i));
			h_Eta_data[i] = (TH1D *)data_file->Get(Form("Eta_AcoOff_%i", i));
			h_Eta_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_Eta_acooff_%i", i));
		}

		if (opt == 6)
		{

			h_FA_mcsignal[i] = (TH1D *)mc_signal->Get(Form("FA_mass_range_%i", i));
			h_FA_mcW[i] = (TH1D *)mc_w->Get(Form("FA_mass_range_%i", i));
			h_FA_mctt[i] = (TH1D *)mc_tt->Get(Form("FA_mass_range_%i", i));
			h_FA_mctau[i] = (TH1D *)mc_signal->Get(Form("FA_tau_mass_range_%i", i));
			h_FA_data[i] = (TH1D *)data_file->Get(Form("FA_mass_range_%i", i));
			h_FA_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_FA_mass_range_%i", i));

			h_Eta_mcsignal[i] = (TH1D *)mc_signal->Get(Form("Eta_mass_range_%i", i));
			h_Eta_mcW[i] = (TH1D *)mc_w->Get(Form("Eta_mass_range_%i", i));
			h_Eta_mctt[i] = (TH1D *)mc_tt->Get(Form("Eta_mass_range_%i", i));
			h_Eta_mctau[i] = (TH1D *)mc_signal->Get(Form("Eta_tau_mass_range_%i", i));
			h_Eta_data[i] = (TH1D *)data_file->Get(Form("Eta_mass_range_%i", i));
			h_Eta_samesign[i] = (TH1D *)mc_signal->Get(Form("samesign_Eta_mass_range_%i", i));
		}
	}

	TH1D *h_DY_weight;
	TH1D *h_tt_weight;
	TH1D *h_W_weight;

	h_DY_weight = (TH1D *)mc_signal->Get("event_weight");
	h_tt_weight = (TH1D *)mc_tt->Get("event_weight");
	h_W_weight = (TH1D *)mc_w->Get("event_weight");

	double w_norm = h_W_weight->Integral();
	double tt_norm = h_tt_weight->Integral();
	double dy_norm = h_DY_weight->Integral();

	plotting_helper *ovo = new plotting_helper();
	ovo->setTDRStyle();

	TH1D *h_normalized_data_FA[11];
	TH1D *h_normalized_mc_FA[11];
	TH1D *h_normalized_mc_bk_FA[11];
	TH1D *h_total_mc_30bins_FA[11];
	TH1D *h_total_mc_120bins_FA[11];

	TH1D *h_normalized_data_Eta[11];
	TH1D *h_normalized_mc_Eta[11];
	TH1D *h_normalized_mc_bk_Eta[11];
	TH1D *h_total_mc_30bins_Eta[11];
	TH1D *h_total_mc_120bins_Eta[11];

	for (int i = 0; i < 11; i++)
	{
		ovo->fixnegativebin(h_FA_mcsignal[i]);
		ovo->fixnegativebin(h_FA_mcW[i]);
		ovo->fixnegativebin(h_FA_mctt[i]);
		ovo->fixnegativebin(h_FA_mctau[i]);
		ovo->fixnegativebin(h_FA_data[i]);
		ovo->fixnegativebin(h_FA_samesign[i]);

		ovo->fixnegativebin(h_Eta_mcsignal[i]);
		ovo->fixnegativebin(h_Eta_mcW[i]);
		ovo->fixnegativebin(h_Eta_mctt[i]);
		ovo->fixnegativebin(h_Eta_mctau[i]);
		ovo->fixnegativebin(h_Eta_data[i]);
		ovo->fixnegativebin(h_Eta_samesign[i]);

		ovo->luminormalize(h_FA_mcW[i], 2, w_norm);
		ovo->luminormalize(h_FA_mctt[i], 3, tt_norm);
		ovo->luminormalize(h_FA_mctau[i], 1, dy_norm);
		ovo->luminormalize(h_FA_mcsignal[i], 1, dy_norm);
		ovo->luminormalize(h_FA_samesign[i], 1, dy_norm);

		ovo->luminormalize(h_Eta_mcW[i], 2, w_norm);
		ovo->luminormalize(h_Eta_mctt[i], 3, tt_norm);
		ovo->luminormalize(h_Eta_mctau[i], 1, dy_norm);
		ovo->luminormalize(h_Eta_mcsignal[i], 1, dy_norm);
		ovo->luminormalize(h_Eta_samesign[i], 1, dy_norm);

		h_normalized_data_FA[i] = (TH1D *)h_FA_data[i]->Clone();
		h_normalized_data_Eta[i] = (TH1D *)h_Eta_data[i]->Clone();

		h_normalized_mc_FA[i] = (TH1D *)h_FA_mcsignal[i]->Clone();
		h_normalized_mc_FA[i]->Add(h_FA_mcW[i], 1);
		h_normalized_mc_FA[i]->Add(h_FA_mctt[i], 1);
		h_normalized_mc_FA[i]->Add(h_FA_mctau[i], 1);
		h_normalized_mc_FA[i]->Add(h_FA_samesign[i], 1);

		h_normalized_mc_Eta[i] = (TH1D *)h_Eta_mcsignal[i]->Clone();
		h_normalized_mc_Eta[i]->Add(h_Eta_mcW[i], 1);
		h_normalized_mc_Eta[i]->Add(h_Eta_mctt[i], 1);
		h_normalized_mc_Eta[i]->Add(h_Eta_mctau[i], 1);
		h_normalized_mc_Eta[i]->Add(h_Eta_samesign[i], 1);

		h_normalized_mc_bk_FA[i] = (TH1D *)h_FA_mcW[i]->Clone();
		h_normalized_mc_bk_FA[i]->Add(h_FA_mctt[i], 1);
		h_normalized_mc_bk_FA[i]->Add(h_FA_mctau[i], 1);
		h_normalized_mc_bk_FA[i]->Add(h_FA_samesign[i], 1);

		h_normalized_mc_bk_Eta[i] = (TH1D *)h_Eta_mcW[i]->Clone();
		h_normalized_mc_bk_Eta[i]->Add(h_Eta_mctt[i], 1);
		h_normalized_mc_bk_Eta[i]->Add(h_Eta_mctau[i], 1);
		h_normalized_mc_bk_Eta[i]->Add(h_Eta_samesign[i], 1);

		h_total_mc_120bins_FA[i] = (TH1D *)h_FA_mcW[i]->Clone();
		h_total_mc_120bins_FA[i]->Add(h_FA_mctt[i]);
		h_total_mc_120bins_FA[i]->Add(h_FA_mctau[i]);
		h_total_mc_120bins_FA[i]->Add(h_FA_mcsignal[i]);
		h_total_mc_120bins_FA[i]->Add(h_FA_samesign[i]);

		h_total_mc_120bins_Eta[i] = (TH1D *)h_Eta_mcW[i]->Clone();
		h_total_mc_120bins_Eta[i]->Add(h_Eta_mctt[i]);
		h_total_mc_120bins_Eta[i]->Add(h_Eta_mctau[i]);
		h_total_mc_120bins_Eta[i]->Add(h_Eta_mcsignal[i]);
		h_total_mc_120bins_Eta[i]->Add(h_Eta_samesign[i]);

		ovo->areanormalize(h_normalized_data_FA[i]);
		ovo->areanormalize(h_normalized_data_Eta[i]);
		ovo->areanormalize(h_normalized_mc_FA[i]);
		ovo->areanormalize(h_normalized_mc_Eta[i]);
		// ovo->areanormalize(h_normalized_mc_bk[i]);

		double normalization_factor_mc_120bins_FA = h_total_mc_120bins_FA[i]->Integral("width");
		double normalization_factor_mc_120bins_Eta = h_total_mc_120bins_Eta[i]->Integral("width");
		h_normalized_mc_bk_FA[i]->Scale(1 / normalization_factor_mc_120bins_FA);
		h_normalized_mc_bk_Eta[i]->Scale(1 / normalization_factor_mc_120bins_Eta);

		h_FA_mcW[i]->Rebin(4);
		h_FA_mctt[i]->Rebin(4);
		h_FA_mctau[i]->Rebin(4);
		h_FA_mcsignal[i]->Rebin(4);
		h_FA_data[i]->Rebin(4);
		h_FA_samesign[i]->Rebin(4);

		h_Eta_mcW[i]->Rebin(4);
		h_Eta_mctt[i]->Rebin(4);
		h_Eta_mctau[i]->Rebin(4);
		h_Eta_mcsignal[i]->Rebin(4);
		h_Eta_data[i]->Rebin(4);
		h_Eta_samesign[i]->Rebin(4);

		h_total_mc_30bins_FA[i] = (TH1D *)h_FA_mcW[i]->Clone();
		h_total_mc_30bins_FA[i]->Add(h_FA_mctt[i]);
		h_total_mc_30bins_FA[i]->Add(h_FA_mctau[i]);
		h_total_mc_30bins_FA[i]->Add(h_FA_mcsignal[i]);
		h_total_mc_30bins_FA[i]->Add(h_FA_samesign[i]);

		h_total_mc_30bins_Eta[i] = (TH1D *)h_Eta_mcW[i]->Clone();
		h_total_mc_30bins_Eta[i]->Add(h_Eta_mctt[i]);
		h_total_mc_30bins_Eta[i]->Add(h_Eta_mctau[i]);
		h_total_mc_30bins_Eta[i]->Add(h_Eta_mcsignal[i]);
		h_total_mc_30bins_Eta[i]->Add(h_Eta_samesign[i]);

		double normalization_factor_mc_30bins_FA = h_total_mc_30bins_FA[i]->Integral("width");
		double normalization_factor_mc_30bins_Eta = h_total_mc_30bins_Eta[i]->Integral("width");

		h_FA_mcW[i]->Scale(1 / normalization_factor_mc_30bins_FA);
		h_FA_mctt[i]->Scale(1 / normalization_factor_mc_30bins_FA);
		h_FA_mctau[i]->Scale(1 / normalization_factor_mc_30bins_FA);
		h_FA_mcsignal[i]->Scale(1 / normalization_factor_mc_30bins_FA);
		h_FA_samesign[i]->Scale(1 / normalization_factor_mc_30bins_FA);

		h_Eta_mcW[i]->Scale(1 / normalization_factor_mc_30bins_Eta);
		h_Eta_mctt[i]->Scale(1 / normalization_factor_mc_30bins_Eta);
		h_Eta_mctau[i]->Scale(1 / normalization_factor_mc_30bins_Eta);
		h_Eta_mcsignal[i]->Scale(1 / normalization_factor_mc_30bins_Eta);
		h_Eta_samesign[i]->Scale(1 / normalization_factor_mc_30bins_Eta);

		ovo->areanormalize(h_FA_data[i]);
		ovo->areanormalize(h_Eta_data[i]);

		ovo->compositeplot(h_FA_data[i], h_FA_mcsignal[i], h_FA_samesign[i], h_FA_mctau[i], h_FA_mcW[i], h_FA_mctt[i], i, opt, false);
		ovo->compositeplot(h_Eta_data[i], h_Eta_mcsignal[i], h_Eta_samesign[i], h_Eta_mctau[i], h_Eta_mcW[i], h_Eta_mctt[i], i, opt, true);

		if (opt == 1)
		{
			ovo->savehistogram(h_normalized_mc_FA[i], h_normalized_data_FA[i], h_normalized_mc_bk_FA[i], i, f_FA_nominal);
			ovo->savehistogram(h_normalized_mc_Eta[i], h_normalized_data_Eta[i], h_normalized_mc_bk_Eta[i], i, f_Eta_nominal);
		}
		if (opt == 2)
		{
			ovo->savehistogram(h_normalized_mc_FA[i], h_normalized_data_FA[i], h_normalized_mc_bk_FA[i], i, f_FA_tnpU);
			ovo->savehistogram(h_normalized_mc_Eta[i], h_normalized_data_Eta[i], h_normalized_mc_bk_Eta[i], i, f_Eta_tnpU);
		}
		if (opt == 3)
		{
			ovo->savehistogram(h_normalized_mc_FA[i], h_normalized_data_FA[i], h_normalized_mc_bk_FA[i], i, f_FA_tnpD);
			ovo->savehistogram(h_normalized_mc_Eta[i], h_normalized_data_Eta[i], h_normalized_mc_bk_Eta[i], i, f_Eta_tnpD);
		}
		if (opt == 4)
		{
			ovo->savehistogram(h_normalized_mc_FA[i], h_normalized_data_FA[i], h_normalized_mc_bk_FA[i], i, f_FA_Acooff);
			ovo->savehistogram(h_normalized_mc_Eta[i], h_normalized_data_Eta[i], h_normalized_mc_bk_Eta[i], i, f_Eta_Acooff);
		}
		if (opt == 6)
		{
			ovo->savehistogram(h_normalized_mc_FA[i], h_normalized_data_FA[i], h_normalized_mc_bk_FA[i], i, f_FA_mass_range);
			ovo->savehistogram(h_normalized_mc_Eta[i], h_normalized_data_Eta[i], h_normalized_mc_bk_Eta[i], i, f_Eta_mass_range);
		}
	}
}