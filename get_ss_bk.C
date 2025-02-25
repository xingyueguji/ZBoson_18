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

// C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

void get_ss_bk()
{
	TH1::SetDefaultSumw2();
	int nbins = 120;
	TFile *mc_file;
	TFile *data_file;

	mc_file = new TFile("./rootfile/mc_signal.root", "UPDATE");
	data_file = new TFile("./rootfile/data_file.root", "READ");

	// Data
	TH1D *FA_nominal[11];
	TH1D *Eta_nominal[11];

	TH1D *FA_AcoOff[11];
	TH1D *Eta_AcoOff[11];

	TH1D *FA_tnpU[11];
	TH1D *Eta_tnpU[11];

	TH1D *FA_tnpD[11];
	TH1D *Eta_tnpD[11];

	TH1D *FA_mass_range[11];
	TH1D *Eta_mass_range[11];

	TH1D *FA_HF_up[11];
	TH1D *Eta_HF_up[11];

	TH1D *FA_HF_down[11];
	TH1D *Eta_HF_down[11];

	// Data same sign

	TH1D *FA_ss_nominal[11];
	TH1D *Eta_ss_nominal[11];

	TH1D *FA_ss_AcoOff[11];
	TH1D *Eta_ss_AcoOff[11];

	TH1D *FA_ss_tnpU[11];
	TH1D *Eta_ss_tnpU[11];

	TH1D *FA_ss_tnpD[11];
	TH1D *Eta_ss_tnpD[11];

	TH1D *FA_ss_mass_range[11];
	TH1D *Eta_ss_mass_range[11];

	TH1D *FA_ss_HF_up[11];
	TH1D *Eta_ss_HF_up[11];

	TH1D *FA_ss_HF_down[11];
	TH1D *Eta_ss_HF_down[11];

	// MC signal
	TH1D *FA_mc_nominal[11];
	TH1D *Eta_mc_nominal[11];

	TH1D *FA_mc_AcoOff[11];
	TH1D *Eta_mc_AcoOff[11];

	TH1D *FA_mc_mass_range[11];
	TH1D *Eta_mc_mass_range[11];

	// Ratio
	TH1D *h_ratio_FA_nominal[11];
	TH1D *h_ratio_FA_tnpU[11];
	TH1D *h_ratio_FA_tnpD[11];
	TH1D *h_ratio_FA_acooff[11];
	TH1D *h_ratio_FA_mass_range[11];
	TH1D *h_ratio_FA_HF_up[11];
	TH1D *h_ratio_FA_HF_down[11];

	TH1D *h_ratio_Eta_nominal[11];
	TH1D *h_ratio_Eta_tnpU[11];
	TH1D *h_ratio_Eta_tnpD[11];
	TH1D *h_ratio_Eta_acooff[11];
	TH1D *h_ratio_Eta_mass_range[11];
	TH1D *h_ratio_Eta_HF_up[11];
	TH1D *h_ratio_Eta_HF_down[11];

	for (int i = 0; i < 11; i++)
	{
		FA_mc_nominal[i] = (TH1D *)mc_file->Get(Form("FA_nominal_%i", i));
		Eta_mc_nominal[i] = (TH1D *)mc_file->Get(Form("Eta_nominal_%i", i));
		FA_mc_AcoOff[i] = (TH1D *)mc_file->Get(Form("FA_AcoOff_%i", i));
		Eta_mc_AcoOff[i] = (TH1D *)mc_file->Get(Form("Eta_AcoOff_%i", i));
		FA_mc_mass_range[i] = (TH1D *)mc_file->Get(Form("FA_mass_range_%i", i));
		Eta_mc_mass_range[i] = (TH1D *)mc_file->Get(Form("Eta_mass_range_%i", i));

		FA_nominal[i] = (TH1D *)data_file->Get(Form("FA_nominal_%i", i));
		Eta_nominal[i] = (TH1D *)data_file->Get(Form("Eta_nominal_%i", i));
		FA_AcoOff[i] = (TH1D *)data_file->Get(Form("FA_AcoOff_%i", i));
		Eta_AcoOff[i] = (TH1D *)data_file->Get(Form("Eta_AcoOff_%i", i));
		FA_tnpU[i] = (TH1D *)data_file->Get(Form("FA_tnpU_%i", i));
		Eta_tnpU[i] = (TH1D *)data_file->Get(Form("Eta_tnpU_%i", i));
		FA_tnpD[i] = (TH1D *)data_file->Get(Form("FA_tnpD_%i", i));
		Eta_tnpD[i] = (TH1D *)data_file->Get(Form("Eta_tnpD_%i", i));
		FA_mass_range[i] = (TH1D *)data_file->Get(Form("FA_mass_range_%i", i));
		Eta_mass_range[i] = (TH1D *)data_file->Get(Form("Eta_mass_range_%i", i));
		FA_HF_up[i] = (TH1D *)data_file->Get(Form("FA_HF_up_%i", i));
		Eta_HF_up[i] = (TH1D *)data_file->Get(Form("Eta_HF_up_%i", i));
		FA_HF_down[i] = (TH1D *)data_file->Get(Form("FA_HF_down_%i", i));
		Eta_HF_down[i] = (TH1D *)data_file->Get(Form("Eta_HF_down_%i", i));

		FA_ss_nominal[i] = (TH1D *)data_file->Get(Form("FA_ss_nominal_%i", i));
		Eta_ss_nominal[i] = (TH1D *)data_file->Get(Form("Eta_ss_nominal_%i", i));
		FA_ss_AcoOff[i] = (TH1D *)data_file->Get(Form("FA_ss_AcoOff_%i", i));
		Eta_ss_AcoOff[i] = (TH1D *)data_file->Get(Form("Eta_ss_AcoOff_%i", i));
		FA_ss_tnpU[i] = (TH1D *)data_file->Get(Form("FA_ss_tnpU_%i", i));
		Eta_ss_tnpU[i] = (TH1D *)data_file->Get(Form("Eta_ss_tnpU_%i", i));
		FA_ss_tnpD[i] = (TH1D *)data_file->Get(Form("FA_ss_tnpD_%i", i));
		Eta_ss_tnpD[i] = (TH1D *)data_file->Get(Form("Eta_ss_tnpD_%i", i));
		FA_ss_mass_range[i] = (TH1D *)data_file->Get(Form("FA_ss_mass_range_%i", i));
		Eta_ss_mass_range[i] = (TH1D *)data_file->Get(Form("Eta_ss_mass_range_%i", i));
		FA_ss_HF_up[i] = (TH1D *)data_file->Get(Form("FA_ss_HF_up_%i", i));
		Eta_ss_HF_up[i] = (TH1D *)data_file->Get(Form("Eta_ss_HF_up_%i", i));
		FA_ss_HF_down[i] = (TH1D *)data_file->Get(Form("FA_ss_HF_down_%i", i));
		Eta_ss_HF_down[i] = (TH1D *)data_file->Get(Form("Eta_ss_HF_down_%i", i));

		h_ratio_FA_nominal[i] = new TH1D(Form("ratio_FA_nominal_%i", i), "", 120, 60, 120);
		h_ratio_FA_tnpU[i] = new TH1D(Form("ratio_FA_tnpU_%i", i), "", 120, 60, 120);
		h_ratio_FA_tnpD[i] = new TH1D(Form("ratio_FA_tnpD_%i", i), "", 120, 60, 120);
		h_ratio_FA_acooff[i] = new TH1D(Form("ratio_FA_acooff_%i", i), "", 120, 60, 120);
		h_ratio_FA_mass_range[i] = new TH1D(Form("ratio_FA_mass_range_%i", i), "", 80, 70, 110);
		h_ratio_FA_HF_up[i] = new TH1D(Form("ratio_FA_HF_up_%i", i), "", 120, 60, 120);
		h_ratio_FA_HF_down[i] = new TH1D(Form("ratio_FA_HF_down_%i", i), "", 120, 60, 120);

		h_ratio_Eta_nominal[i] = new TH1D(Form("ratio_Eta_nominal_%i", i), "", 120, 60, 120);
		h_ratio_Eta_tnpU[i] = new TH1D(Form("ratio_Eta_tnpU_%i", i), "", 120, 60, 120);
		h_ratio_Eta_tnpD[i] = new TH1D(Form("ratio_Eta_tnpD_%i", i), "", 120, 60, 120);
		h_ratio_Eta_acooff[i] = new TH1D(Form("ratio_Eta_acooff_%i", i), "", 120, 60, 120);
		h_ratio_Eta_mass_range[i] = new TH1D(Form("ratio_Eta_mass_range_%i", i), "", 80, 70, 110);
		h_ratio_Eta_HF_up[i] = new TH1D(Form("ratio_Eta_HF_up_%i", i), "", 120, 60, 120);
		h_ratio_Eta_HF_down[i] = new TH1D(Form("ratio_Eta_HF_down_%i", i), "", 120, 60, 120);
	}

	for (int i = 0; i < 11; i++)
	{
		for (int j = 1; j <= 80; j++)
		{
			double content_mc_os_FA_mass_range = FA_mc_mass_range[i]->GetBinContent(j);
			double content_data_os_FA_mass_range = FA_mass_range[i]->GetBinContent(j);
			double content_data_ss_FA_mass_range = FA_ss_mass_range[i]->GetBinContent(j);
			double ratio_FA_mass_range;

			double content_mc_os_Eta_mass_range = Eta_mc_mass_range[i]->GetBinContent(j);
			double content_data_os_Eta_mass_range = Eta_mass_range[i]->GetBinContent(j);
			double content_data_ss_Eta_mass_range = Eta_ss_mass_range[i]->GetBinContent(j);
			double ratio_Eta_mass_range;

			if (content_data_os_FA_mass_range > 0)
			{
				ratio_FA_mass_range = content_data_ss_FA_mass_range / content_data_os_FA_mass_range;
			}
			else
			{
				ratio_FA_mass_range = 0;
			}

			if (content_data_os_Eta_mass_range > 0)
			{
				ratio_Eta_mass_range = content_data_ss_Eta_mass_range / content_data_os_Eta_mass_range;
			}
			else
			{
				ratio_Eta_mass_range = 0;
			}

			double estimate_mc_ss_FA_mass_range = content_mc_os_FA_mass_range * ratio_FA_mass_range;
			double estimate_mc_ss_Eta_mass_range = content_mc_os_Eta_mass_range * ratio_Eta_mass_range;

			if (estimate_mc_ss_FA_mass_range >= 0)
				h_ratio_FA_mass_range[i]->SetBinContent(j, estimate_mc_ss_FA_mass_range);
			if (estimate_mc_ss_Eta_mass_range >= 0)
				h_ratio_Eta_mass_range[i]->SetBinContent(j, estimate_mc_ss_Eta_mass_range);
		}
		for (int j = 1; j <= 120; j++)
		{

			double content_mc_os_FA_nominal = FA_mc_nominal[i]->GetBinContent(j);
			double content_data_os_FA_nominal = FA_nominal[i]->GetBinContent(j);
			double content_data_ss_FA_nominal = FA_ss_nominal[i]->GetBinContent(j);
			double ratio_FA_nominal;

			// double content_mc_os_FA_tnpU no such thing, used FA_mc_nominal
			double content_data_os_FA_tnpU = FA_tnpU[i]->GetBinContent(j);
			double content_data_ss_FA_tnpU = FA_ss_tnpU[i]->GetBinContent(j);
			double ratio_FA_tnpU;

			// double content_mc_os_FA_tnpD no such thing, used FA_mc_nominal
			double content_data_os_FA_tnpD = FA_tnpD[i]->GetBinContent(j);
			double content_data_ss_FA_tnpD = FA_ss_tnpD[i]->GetBinContent(j);
			double ratio_FA_tnpD;

			double content_mc_os_FA_acooff = FA_mc_AcoOff[i]->GetBinContent(j);
			double content_data_os_FA_acooff = FA_AcoOff[i]->GetBinContent(j);
			double content_data_ss_FA_acooff = FA_ss_AcoOff[i]->GetBinContent(j);
			double ratio_FA_acooff;

			double content_data_os_FA_HF_up = FA_HF_up[i]->GetBinContent(j);
			double content_data_ss_FA_HF_up = FA_ss_HF_up[i]->GetBinContent(j);
			double ratio_FA_HF_up;

			double content_data_os_FA_HF_down = FA_HF_down[i]->GetBinContent(j);
			double content_data_ss_FA_HF_down = FA_ss_HF_down[i]->GetBinContent(j);
			double ratio_FA_HF_down;

			if (content_data_os_FA_nominal > 0)
			{
				ratio_FA_nominal = content_data_ss_FA_nominal / content_data_os_FA_nominal;
			}
			else
			{
				ratio_FA_nominal = 0;
			}

			if (content_data_os_FA_tnpU > 0)
			{
				ratio_FA_tnpU = content_data_ss_FA_tnpU / content_data_os_FA_tnpU;
			}
			else
			{
				ratio_FA_tnpU = 0;
			}

			if (content_data_os_FA_tnpD > 0)
			{
				ratio_FA_tnpD = content_data_ss_FA_tnpD / content_data_os_FA_tnpD;
			}
			else
			{
				ratio_FA_tnpD = 0;
			}

			if (content_data_os_FA_acooff > 0)
			{
				ratio_FA_acooff = content_data_ss_FA_acooff / content_data_os_FA_acooff;
			}
			else
			{
				ratio_FA_acooff = 0;
			}

			if (content_data_os_FA_HF_up > 0)
			{
				ratio_FA_HF_up = content_data_ss_FA_HF_up / content_data_os_FA_HF_up;
			}
			else
			{
				ratio_FA_HF_up = 0;
			}

			if (content_data_os_FA_HF_down > 0)
			{
				ratio_FA_HF_down = content_data_ss_FA_HF_down / content_data_os_FA_HF_down;
			}
			else
			{
				ratio_FA_HF_down = 0;
			}

			double content_mc_os_Eta_nominal = Eta_mc_nominal[i]->GetBinContent(j);
			double content_data_os_Eta_nominal = Eta_nominal[i]->GetBinContent(j);
			double content_data_ss_Eta_nominal = Eta_ss_nominal[i]->GetBinContent(j);
			double ratio_Eta_nominal;

			// double content_mc_os_Eta_tnpU no such thing, used Eta_mc_nominal
			double content_data_os_Eta_tnpU = Eta_tnpU[i]->GetBinContent(j);
			double content_data_ss_Eta_tnpU = Eta_ss_tnpU[i]->GetBinContent(j);
			double ratio_Eta_tnpU;

			// double content_mc_os_Eta_tnpD no such thing, used Eta_mc_nominal
			double content_data_os_Eta_tnpD = Eta_tnpD[i]->GetBinContent(j);
			double content_data_ss_Eta_tnpD = Eta_ss_tnpD[i]->GetBinContent(j);
			double ratio_Eta_tnpD;

			double content_mc_os_Eta_acooff = Eta_mc_AcoOff[i]->GetBinContent(j);
			double content_data_os_Eta_acooff = Eta_AcoOff[i]->GetBinContent(j);
			double content_data_ss_Eta_acooff = Eta_ss_AcoOff[i]->GetBinContent(j);
			double ratio_Eta_acooff;

			double content_data_os_Eta_HF_up = Eta_HF_up[i]->GetBinContent(j);
			double content_data_ss_Eta_HF_up = Eta_ss_HF_up[i]->GetBinContent(j);
			double ratio_Eta_HF_up;

			double content_data_os_Eta_HF_down = Eta_HF_down[i]->GetBinContent(j);
			double content_data_ss_Eta_HF_down = Eta_ss_HF_down[i]->GetBinContent(j);
			double ratio_Eta_HF_down;

			if (content_data_os_Eta_nominal > 0)
			{
				ratio_Eta_nominal = content_data_ss_Eta_nominal / content_data_os_Eta_nominal;
			}
			else
			{
				ratio_Eta_nominal = 0;
			}

			if (content_data_os_Eta_tnpU > 0)
			{
				ratio_Eta_tnpU = content_data_ss_Eta_tnpU / content_data_os_Eta_tnpU;
			}
			else
			{
				ratio_Eta_tnpU = 0;
			}

			if (content_data_os_Eta_tnpD > 0)
			{
				ratio_Eta_tnpD = content_data_ss_Eta_tnpD / content_data_os_Eta_tnpD;
			}
			else
			{
				ratio_Eta_tnpD = 0;
			}

			if (content_data_os_Eta_acooff > 0)
			{
				ratio_Eta_acooff = content_data_ss_Eta_acooff / content_data_os_Eta_acooff;
			}
			else
			{
				ratio_Eta_acooff = 0;
			}

			if (content_data_os_Eta_HF_up > 0)
			{
				ratio_Eta_HF_up = content_data_ss_Eta_HF_up / content_data_os_Eta_HF_up;
			}
			else
			{
				ratio_Eta_HF_up = 0;
			}

			if (content_data_os_Eta_HF_down > 0)
			{
				ratio_Eta_HF_down = content_data_ss_Eta_HF_down / content_data_os_Eta_HF_down;
			}
			else
			{
				ratio_Eta_HF_down = 0;
			}

			double estimate_mc_ss_FA_nominal = content_mc_os_FA_nominal * ratio_FA_nominal;
			double estimate_mc_ss_FA_tnpU = content_mc_os_FA_nominal * ratio_FA_tnpU;
			double estimate_mc_ss_FA_tnpD = content_mc_os_FA_nominal * ratio_FA_tnpD;
			double estimate_mc_ss_FA_acooff = content_mc_os_FA_acooff * ratio_FA_acooff;
			double estimate_mc_ss_FA_HF_up = content_mc_os_FA_nominal * ratio_FA_HF_up;
			double estimate_mc_ss_FA_HF_down = content_mc_os_FA_nominal * ratio_FA_HF_down;

			double estimate_mc_ss_Eta_nominal = content_mc_os_Eta_nominal * ratio_Eta_nominal;
			double estimate_mc_ss_Eta_tnpU = content_mc_os_Eta_nominal * ratio_Eta_tnpU;
			double estimate_mc_ss_Eta_tnpD = content_mc_os_Eta_nominal * ratio_Eta_tnpD;
			double estimate_mc_ss_Eta_acooff = content_mc_os_Eta_acooff * ratio_Eta_acooff;
			double estimate_mc_ss_Eta_HF_up = content_mc_os_Eta_nominal * ratio_Eta_HF_up;
			double estimate_mc_ss_Eta_HF_down = content_mc_os_Eta_nominal * ratio_Eta_HF_down;

			if (estimate_mc_ss_FA_nominal >= 0)
				h_ratio_FA_nominal[i]->SetBinContent(j, estimate_mc_ss_FA_nominal);
			if (estimate_mc_ss_FA_tnpU >= 0)
				h_ratio_FA_tnpU[i]->SetBinContent(j, estimate_mc_ss_FA_tnpU);
			if (estimate_mc_ss_FA_tnpD >= 0)
				h_ratio_FA_tnpD[i]->SetBinContent(j, estimate_mc_ss_FA_tnpD);
			if (estimate_mc_ss_FA_acooff >= 0)
				h_ratio_FA_acooff[i]->SetBinContent(j, estimate_mc_ss_FA_acooff);
			if (estimate_mc_ss_FA_HF_up >= 0)
				h_ratio_FA_HF_up[i]->SetBinContent(j, estimate_mc_ss_FA_HF_up);
			if (estimate_mc_ss_FA_HF_down >= 0)
				h_ratio_FA_HF_down[i]->SetBinContent(j, estimate_mc_ss_FA_HF_down);

			if (estimate_mc_ss_Eta_nominal >= 0)
				h_ratio_Eta_nominal[i]->SetBinContent(j, estimate_mc_ss_Eta_nominal);
			if (estimate_mc_ss_Eta_tnpU >= 0)
				h_ratio_Eta_tnpU[i]->SetBinContent(j, estimate_mc_ss_Eta_tnpU);
			if (estimate_mc_ss_Eta_tnpD >= 0)
				h_ratio_Eta_tnpD[i]->SetBinContent(j, estimate_mc_ss_Eta_tnpD);
			if (estimate_mc_ss_Eta_acooff >= 0)
				h_ratio_Eta_acooff[i]->SetBinContent(j, estimate_mc_ss_Eta_acooff);
			if (estimate_mc_ss_Eta_HF_up >= 0)
				h_ratio_Eta_HF_up[i]->SetBinContent(j, estimate_mc_ss_Eta_HF_up);
			if (estimate_mc_ss_Eta_HF_down >= 0)
				h_ratio_Eta_HF_down[i]->SetBinContent(j, estimate_mc_ss_Eta_HF_down);

			// if (i == 7) cout << "Bin content is " <<  estimate_mc_ss << endl;
		}
	}

	mc_file->cd();
	for (int i = 0; i < 11; i++)
	{
		h_ratio_FA_nominal[i]->Write(Form("samesign_FA_nominal_%i", i), 2);
		h_ratio_FA_tnpU[i]->Write(Form("samesign_FA_tnpU_%i", i), 2);
		h_ratio_FA_tnpD[i]->Write(Form("samesign_FA_tnpD_%i", i), 2);
		h_ratio_FA_acooff[i]->Write(Form("samesign_FA_acooff_%i", i), 2);
		h_ratio_FA_HF_up[i]->Write(Form("samesign_FA_HF_up_%i", i), 2);
		h_ratio_FA_HF_down[i]->Write(Form("samesign_FA_HF_down_%i", i), 2);

		h_ratio_Eta_nominal[i]->Write(Form("samesign_Eta_nominal_%i", i), 2);
		h_ratio_Eta_tnpU[i]->Write(Form("samesign_Eta_tnpU_%i", i), 2);
		h_ratio_Eta_tnpD[i]->Write(Form("samesign_Eta_tnpD_%i", i), 2);
		h_ratio_Eta_acooff[i]->Write(Form("samesign_Eta_acooff_%i", i), 2);
		h_ratio_Eta_HF_up[i]->Write(Form("samesign_Eta_HF_up_%i", i), 2);
		h_ratio_Eta_HF_down[i]->Write(Form("samesign_Eta_HF_down_%i", i), 2);

		h_ratio_FA_mass_range[i]->Write(Form("samesign_FA_mass_range_%i", i), 2);
		h_ratio_Eta_mass_range[i]->Write(Form("samesign_Eta_mass_range_%i", i), 2);
	}
	mc_file->Close();
	data_file->Close();
}