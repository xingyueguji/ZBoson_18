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

void get_data()
{
	TH1::SetDefaultSumw2();
	gSystem->Load("./header/libDict.so");

	MC_18 *data = new MC_18();
	MC_18 *data_same_sign = new MC_18();

	data->SetupRootfile(4, 0);
	data->SetupBranches(1);
	data_same_sign->SetupRootfile(5, 0);
	data_same_sign->SetupBranches(1);

	MCReweight *vzRW = new MCReweight("WeightsAndEfficiencies_Z2mummu/vzReweight.root", "WeightsAndEfficiencies_Z2mummu/centralityFlatteningWeight.root");
	PtReweightSpectrum *spectrumRW = new PtReweightSpectrum("WeightsAndEfficiencies_Z2mummu/ptSpectrumReweighting.root");
	MuonTnP *tnp = new MuonTnP();

	TEfficiency *e[11];
	TEfficiency *e_up[11];
	TEfficiency *e_down[11];
	TEfficiency *e_acooff[11];

	TFile *eff_f1 = new TFile("./rootfile/mc_eff.root", "READ");
	for (int i = 0; i < 11; i++)
	{
		e[i] = (TEfficiency *)eff_f1->Get(Form("eff_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]));
		e_up[i] = (TEfficiency *)eff_f1->Get(Form("eff_U_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]));
		e_down[i] = (TEfficiency *)eff_f1->Get(Form("eff_D_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]));
		e_acooff[i] = (TEfficiency *)eff_f1->Get(Form("eff_noAco_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]));
	}

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

	TH1D *FA_HF_up[11];
	TH1D *Eta_HF_up[11];

	TH1D *FA_HF_down[11];
	TH1D *Eta_HF_down[11];

	TH1D *FA_ss_HF_up[11];
	TH1D *Eta_ss_HF_up[11];

	TH1D *FA_ss_HF_down[11];
	TH1D *Eta_ss_HF_down[11];

	TH1D *pT_spec_FA = new TH1D("pT_spec_FA", "", 200, 0, 200);
	TH1D *pT_spec_Eta = new TH1D("pT_spec_Eta", "", 200, 0, 200);

	for (int i = 0; i < data->centarraysize; i++)
	{
		FA_nominal[i] = new TH1D(Form("FA_nominal_%i", i), "", 120, 60, 120);
		Eta_nominal[i] = new TH1D(Form("Eta_nominal_%i", i), "", 120, 60, 120);
		FA_AcoOff[i] = new TH1D(Form("FA_AcoOff_%i", i), "", 120, 60, 120);
		Eta_AcoOff[i] = new TH1D(Form("Eta_AcoOff_%i", i), "", 120, 60, 120);
		FA_tnpU[i] = new TH1D(Form("FA_tnpU_%i", i), "", 120, 60, 120);
		Eta_tnpU[i] = new TH1D(Form("Eta_tnpU_%i", i), "", 120, 60, 120);
		FA_tnpD[i] = new TH1D(Form("FA_tnpD_%i", i), "", 120, 60, 120);
		Eta_tnpD[i] = new TH1D(Form("Eta_tnpD_%i", i), "", 120, 60, 120);
		FA_mass_range[i] = new TH1D(Form("FA_mass_range_%i", i), "", 80, 70, 110);
		Eta_mass_range[i] = new TH1D(Form("Eta_mass_range_%i", i), "", 80, 70, 110);
		FA_HF_up[i] = new TH1D(Form("FA_HF_up_%i", i), "", 120, 60, 120);
		Eta_HF_up[i] = new TH1D(Form("Eta_HF_up_%i", i), "", 120, 60, 120);
		FA_HF_down[i] = new TH1D(Form("FA_HF_down_%i", i), "", 120, 60, 120);
		Eta_HF_down[i] = new TH1D(Form("Eta_HF_down_%i", i), "", 120, 60, 120);

		FA_ss_nominal[i] = new TH1D(Form("FA_ss_nominal_%i", i), "", 120, 60, 120);
		Eta_ss_nominal[i] = new TH1D(Form("Eta_ss_nominal_%i", i), "", 120, 60, 120);
		FA_ss_AcoOff[i] = new TH1D(Form("FA_ss_AcoOff_%i", i), "", 120, 60, 120);
		Eta_ss_AcoOff[i] = new TH1D(Form("Eta_ss_AcoOff_%i", i), "", 120, 60, 120);
		FA_ss_tnpU[i] = new TH1D(Form("FA_ss_tnpU_%i", i), "", 120, 60, 120);
		Eta_ss_tnpU[i] = new TH1D(Form("Eta_ss_tnpU_%i", i), "", 120, 60, 120);
		FA_ss_tnpD[i] = new TH1D(Form("FA_ss_tnpD_%i", i), "", 120, 60, 120);
		Eta_ss_tnpD[i] = new TH1D(Form("Eta_ss_tnpD_%i", i), "", 120, 60, 120);
		FA_ss_mass_range[i] = new TH1D(Form("FA_ss_mass_range_%i", i), "", 80, 70, 110);
		Eta_ss_mass_range[i] = new TH1D(Form("Eta_ss_mass_range_%i", i), "", 80, 70, 110);
		FA_ss_HF_up[i] = new TH1D(Form("FA_ss_HF_up_%i", i), "", 120, 60, 120);
		Eta_ss_HF_up[i] = new TH1D(Form("Eta_ss_HF_up_%i", i), "", 120, 60, 120);
		FA_ss_HF_down[i] = new TH1D(Form("FA_ss_HF_down_%i", i), "", 120, 60, 120);
		Eta_ss_HF_down[i] = new TH1D(Form("Eta_ss_HF_down_%i", i), "", 120, 60, 120);
	}

	int nentries = data->t1->GetEntries();
	cout << "We have " << nentries << " events in total" << endl;

	for (int i = 0; i < nentries; i++)
	{
		if (i % 100000 == 0)
			cout << "Event " << i << endl;
		// cout << "Event " << i << endl;
		data->t1->GetEntry(i);

		// Cut setup
		if (!(data->evtSel[1]))
			continue;
		if (!(data->evtSel[2]))
			continue;
		if (!(data->evtSel[3]))
			continue;
		if (abs(data->bestvtxZ) > 15)
			continue;

		int hiBin = 0;
		int hiBin_UP = 0;
		int hiBin_DOWN = 0;

		hiBin = data->getcentrality(0);
		hiBin_UP = data->getcentrality(1);
		hiBin_DOWN = data->getcentrality(2);

		if (!(data->trigHLT[6]))
			continue;

		for (int j = 0; j < data->candSize; j++)
		{
			if (data->mass[j] < data->masslowlimit || data->mass[j] > data->masshighlimit)
				continue;
			if (data->pTD1[j] < data->ptlowlimit)
				continue;
			if (data->pTD2[j] < data->ptlowlimit)
				continue;
			if (abs(data->EtaD1[j]) > 2.4)
				continue;
			if (abs(data->EtaD2[j]) > 2.4)
				continue;
			if (data->VtxProb[j] < 0.001)
				continue;
			if (abs(data->y[j]) > 2.4)
				continue;
			if (!(data->tightMuon1Cut(j, "POG") && data->tightMuon2Cut(j, "POG")))
				continue;

			bool isDaughter1Trigger = data->trigMuon1->at(6).at(j); // 6 is HLT_HIL3Mu12;
			bool isDaughter2Trigger = data->trigMuon2->at(6).at(j);
			if (!(isDaughter1Trigger || isDaughter2Trigger))
				continue;

			bool isOppositeSign = data->chargeD1[j] != data->chargeD2[j];

			if (!isOppositeSign)
				continue;

			int centbinpositioncounter[11] = {};
			data->CentBinSearching(centbinpositioncounter, hiBin);
			int centbinpositioncounter_UP[11] = {};
			data->CentBinSearching(centbinpositioncounter_UP, hiBin_UP);
			int centbinpositioncounter_DOWN[11] = {};
			data->CentBinSearching(centbinpositioncounter_DOWN, hiBin_DOWN);

			float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(data->PhiD1[j] - data->PhiD2[j]))) / TMath::Pi();
			bool passesAco[3] = {1, 1, 1};
			if (data->pT[j] < 1.25 && acoplanarity < 0.001)
				passesAco[0] = false;

			// Here's for generate pT spectrum of PbPb

			if (passesAco[0])
			{
				double efficiency = data->getEfficiency(e[10], data->y[j], data->pT[j]);
				pT_spec_FA->Fill(data->pT[j], 1.0 / efficiency);
				if ((abs(data->EtaD1[j]) < 1) && (abs(data->EtaD2[j]) < 1))
				{
					double efficiency = data->getEfficiency(e[10], data->y[j], data->pT[j]);
					pT_spec_Eta->Fill(data->pT[j], 1.0 / efficiency);
				}
			}

			for (int k = 0; k < data->centarraysize; k++)
			{
				if (centbinpositioncounter[k] != 0)
				{
					if (passesAco[0])
					{
						double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
						double efficiency_up = data->getEfficiency(e_up[k], data->y[j], data->pT[j]);
						double efficiency_down = data->getEfficiency(e_down[k], data->y[j], data->pT[j]);

						FA_nominal[k]->Fill(data->mass[j], 1.0 / efficiency);
						FA_tnpU[k]->Fill(data->mass[j], 1.0 / efficiency_up);
						FA_tnpD[k]->Fill(data->mass[j], 1.0 / efficiency_down);
						FA_mass_range[k]->Fill(data->mass[j], 1.0 / efficiency);

						// with eta cut < 1 on both D1 and D2 muons
						if ((abs(data->EtaD1[j]) < 1) && (abs(data->EtaD2[j]) < 1))
						{
							double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
							double efficiency_up = data->getEfficiency(e_up[k], data->y[j], data->pT[j]);
							double efficiency_down = data->getEfficiency(e_down[k], data->y[j], data->pT[j]);
							Eta_nominal[k]->Fill(data->mass[j], 1.0 / efficiency);
							Eta_tnpU[k]->Fill(data->mass[j], 1.0 / efficiency_up);
							Eta_tnpD[k]->Fill(data->mass[j], 1.0 / efficiency_down);
							Eta_mass_range[k]->Fill(data->mass[j], 1.0 / efficiency);
						}
					}

					double efficiency_acooff = data->getEfficiency(e_acooff[k], data->y[j], data->pT[j]);
					FA_AcoOff[k]->Fill(data->mass[j], 1.0 / efficiency_acooff);

					if ((abs(data->EtaD1[j]) < 1) && (abs(data->EtaD2[j]) < 1))
					{
						double efficiency_acooff = data->getEfficiency(e_acooff[k], data->y[j], data->pT[j]);
						Eta_AcoOff[k]->Fill(data->mass[j], 1.0 / efficiency_acooff);
					}
				}
				if (centbinpositioncounter_UP[k] != 0)
				{
					if (passesAco[0])
					{
						double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
						FA_HF_up[k]->Fill(data->mass[j], 1.0 / efficiency);

						// with eta cut < 1 on both D1 and D2 muons
						if ((abs(data->EtaD1[j]) < 1) && (abs(data->EtaD2[j]) < 1))
						{
							double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
							Eta_HF_up[k]->Fill(data->mass[j], 1.0 / efficiency);
						}
					}
				}
				if (centbinpositioncounter_DOWN[k] != 0)
				{
					if (passesAco[0])
					{
						double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
						FA_HF_down[k]->Fill(data->mass[j], 1.0 / efficiency);

						// with eta cut < 1 on both D1 and D2 muons
						if ((abs(data->EtaD1[j]) < 1) && (abs(data->EtaD2[j]) < 1))
						{
							double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
							Eta_HF_down[k]->Fill(data->mass[j], 1.0 / efficiency);
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < data_same_sign->t1->GetEntries(); i++)
	{
		if (i % 100000 == 0)
			cout << "Event " << i << endl;
		// cout << "Event " << i << endl;
		data_same_sign->t1->GetEntry(i);

		// Cut setup
		if (!(data_same_sign->evtSel[1]))
			continue;
		if (!(data_same_sign->evtSel[2]))
			continue;
		if (!(data_same_sign->evtSel[3]))
			continue;
		if (abs(data_same_sign->bestvtxZ) > 15)
			continue;

		int hiBin = 0;
		int hiBin_UP = 0;
		int hiBin_DOWN = 0;

		hiBin = data_same_sign->getcentrality(0); // FIXME: doZDC? hiBinVar?
		hiBin_UP = data_same_sign->getcentrality(1);
		hiBin_DOWN = data_same_sign->getcentrality(2);

		int centbinpositioncounter[11] = {};
		data_same_sign->CentBinSearching(centbinpositioncounter, hiBin);
		int centbinpositioncounter_UP[11] = {};
		data_same_sign->CentBinSearching(centbinpositioncounter_UP, hiBin_UP);
		int centbinpositioncounter_DOWN[11] = {};
		data_same_sign->CentBinSearching(centbinpositioncounter_DOWN, hiBin_DOWN);

		if (!(data_same_sign->trigHLT[6]))
			continue;

		for (int j = 0; j < data_same_sign->candSize; j++)
		{
			if (data_same_sign->mass[j] < data_same_sign->masslowlimit || data_same_sign->mass[j] > data_same_sign->masshighlimit)
				continue;
			if (data_same_sign->pTD1[j] < data_same_sign->ptlowlimit)
				continue;
			if (data_same_sign->pTD2[j] < data_same_sign->ptlowlimit)
				continue;
			if (abs(data_same_sign->EtaD1[j]) > 2.4)
				continue;
			if (abs(data_same_sign->EtaD2[j]) > 2.4)
				continue;
			if (data_same_sign->VtxProb[j] < 0.001)
				continue;
			if (abs(data_same_sign->y[j]) > 2.4)
				continue;
			if (!(data_same_sign->tightMuon1Cut(j, "POG") && data_same_sign->tightMuon2Cut(j, "POG")))
				continue;

			bool isDaughter1Trigger = data_same_sign->trigMuon1->at(6).at(j); // 6 is HLT_HIL3Mu12;
			bool isDaughter2Trigger = data_same_sign->trigMuon2->at(6).at(j);
			if (!(isDaughter1Trigger || isDaughter2Trigger))
				continue;

			bool isOppositeSign = data_same_sign->chargeD1[j] != data_same_sign->chargeD2[j];

			if (isOppositeSign)
				continue;

			int centbinpositioncounter[11] = {};
			data_same_sign->CentBinSearching(centbinpositioncounter, hiBin);

			float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(data_same_sign->PhiD1[j] - data_same_sign->PhiD2[j]))) / TMath::Pi();
			bool passesAco[3] = {1, 1, 1};
			if (data_same_sign->pT[j] < 1.25 && acoplanarity < 0.001)
				passesAco[0] = false;

			for (int k = 0; k < data_same_sign->centarraysize; k++)
			{
				if (centbinpositioncounter[k] != 0)
				{
					if (passesAco[0])
					{
						double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
						double efficiency_up = data_same_sign->getEfficiency(e_up[k], data_same_sign->y[j], data_same_sign->pT[j]);
						double efficiency_down = data_same_sign->getEfficiency(e_down[k], data_same_sign->y[j], data_same_sign->pT[j]);

						FA_ss_nominal[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
						FA_ss_tnpU[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency_up);
						FA_ss_tnpD[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency_down);
						FA_ss_mass_range[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);

						// with eta cut < 1 on both D1 and D2 muons
						if ((abs(data_same_sign->EtaD1[j]) < 1) && (abs(data_same_sign->EtaD2[j]) < 1))
						{
							double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
							double efficiency_up = data_same_sign->getEfficiency(e_up[k], data_same_sign->y[j], data_same_sign->pT[j]);
							double efficiency_down = data_same_sign->getEfficiency(e_down[k], data_same_sign->y[j], data_same_sign->pT[j]);
							Eta_ss_nominal[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
							Eta_ss_tnpU[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency_up);
							Eta_ss_tnpD[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency_down);
							Eta_ss_mass_range[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
						}
					}

					double efficiency_acooff = data_same_sign->getEfficiency(e_acooff[k], data_same_sign->y[j], data_same_sign->pT[j]);
					FA_ss_AcoOff[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency_acooff);

					if ((abs(data_same_sign->EtaD1[j]) < 1) && (abs(data_same_sign->EtaD2[j]) < 1))
					{
						double efficiency_acooff = data_same_sign->getEfficiency(e_acooff[k], data_same_sign->y[j], data_same_sign->pT[j]);
						Eta_ss_AcoOff[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency_acooff);
					}
				}
				if (centbinpositioncounter_UP[k] != 0)
				{
					if (passesAco[0])
					{
						double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
						FA_ss_HF_up[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);

						// with eta cut < 1 on both D1 and D2 muons
						if ((abs(data_same_sign->EtaD1[j]) < 1) && (abs(data_same_sign->EtaD2[j]) < 1))
						{
							double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
							Eta_ss_HF_up[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
						}
					}
				}
				if (centbinpositioncounter_DOWN[k] != 0)
				{
					if (passesAco[0])
					{
						double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
						FA_ss_HF_down[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);

						// with eta cut < 1 on both D1 and D2 muons
						if ((abs(data_same_sign->EtaD1[j]) < 1) && (abs(data_same_sign->EtaD2[j]) < 1))
						{
							double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
							Eta_ss_HF_down[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
						}
					}
				}
			}
		}
	}

	TFile *histogram_file;
	histogram_file = new TFile("./rootfile/data_file.root", "UPDATE");
	histogram_file->cd();
	for (int i = 0; i < data->centarraysize; i++)
	{
		FA_nominal[i]->Write("", 2);
		Eta_nominal[i]->Write("", 2);
		FA_AcoOff[i]->Write("", 2);
		Eta_AcoOff[i]->Write("", 2);
		FA_tnpU[i]->Write("", 2);
		Eta_tnpU[i]->Write("", 2);
		FA_tnpD[i]->Write("", 2);
		Eta_tnpD[i]->Write("", 2);
		FA_mass_range[i]->Write("", 2);
		Eta_mass_range[i]->Write("", 2);
		FA_HF_up[i]->Write("", 2);
		Eta_HF_up[i]->Write("", 2);
		FA_HF_down[i]->Write("", 2);
		Eta_HF_down[i]->Write("", 2);

		FA_ss_nominal[i]->Write("", 2);
		Eta_ss_nominal[i]->Write("", 2);
		FA_ss_AcoOff[i]->Write("", 2);
		Eta_ss_AcoOff[i]->Write("", 2);
		FA_ss_tnpU[i]->Write("", 2);
		Eta_ss_tnpU[i]->Write("", 2);
		FA_ss_tnpD[i]->Write("", 2);
		Eta_ss_tnpD[i]->Write("", 2);
		FA_ss_mass_range[i]->Write("", 2);
		Eta_ss_mass_range[i]->Write("", 2);
		FA_ss_HF_up[i]->Write("", 2);
		Eta_ss_HF_up[i]->Write("", 2);
		FA_ss_HF_down[i]->Write("", 2);
		Eta_ss_HF_down[i]->Write("", 2);
	}
	histogram_file->Close();

	TFile *pT_file = new TFile("./rootfile/pT_file.root", "UPDATE");
	pT_file->cd();
	pT_spec_FA->Write("", 2);
	pT_spec_Eta->Write("", 2);

	data->f1->Close();
}