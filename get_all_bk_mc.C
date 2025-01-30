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

void get_all_bk_mc(int opt = 3)
{

	auto start = std::chrono::high_resolution_clock::now();
	TH1::SetDefaultSumw2();
	// opt == 1 means signal MC;
	// opt == 2 means W;
	// opt == 3 means tt;
	// dont have to worry about starlight

	gSystem->Load("./header/libDict.so");

	MC_18 *s = new MC_18();

	MCReweight *vzRW = new MCReweight("WeightsAndEfficiencies_Z2mummu/vzReweight.root", "WeightsAndEfficiencies_Z2mummu/centralityFlatteningWeight.root");
	PtReweightSpectrum *spectrumRW = new PtReweightSpectrum("WeightsAndEfficiencies_Z2mummu/ptSpectrumReweighting.root");
	MuonTnP *tnp = new MuonTnP();

	if (opt == 1)
	{
		s->SetupRootfile(1, 1);
		s->SetupBranches(0);
	}
	if (opt == 2)
	{
		s->SetupRootfile(2, 1);
		s->SetupBranches(0);
	}
	if (opt == 3)
	{
		s->SetupRootfile(3, 1);
		s->SetupBranches(0);
	}

	TEfficiency *e[11];
	TEfficiency *e_without_aco[11];

	TFile *eff_f1 = new TFile("./rootfile/mc_eff.root", "READ");
	for (int i = 0; i < 11; i++)
	{
		e[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
		e_without_aco[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_noAco_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
	}

	// Histogram setup

	// These are nominal, so with Aco cut, with eff.
	TH1D *event_weight = new TH1D("event_weight", "event_weight", 1, 0, 1);

	TH1D *FA_nominal[11];
	TH1D *Eta_nominal[11];

	TH1D *FA_AcoOff[11];
	TH1D *Eta_AcoOff[11];

	TH1D *FA_tau_nominal[11];
	TH1D *Eta_tau_nominal[11];

	TH1D *FA_tau_AcoOff[11];
	TH1D *Eta_tau_AcoOff[11];

	for (int i = 0; i < s->centarraysize; i++)
	{
		FA_nominal[i] = new TH1D(Form("FA_nominal_%i", i), "", 120, 60, 120);
		Eta_nominal[i] = new TH1D(Form("Eta_nominal_%i", i), "", 120, 60, 120);
		FA_AcoOff[i] = new TH1D(Form("FA_AcoOff_%i", i), "", 120, 60, 120);
		Eta_AcoOff[i] = new TH1D(Form("Eta_AcoOff_%i", i), "", 120, 60, 120);

		FA_tau_nominal[i] = new TH1D(Form("FA_tau_nominal_%i", i), "", 120, 60, 120);
		Eta_tau_nominal[i] = new TH1D(Form("Eta_tau_nominal_%i", i), "", 120, 60, 120);
		FA_tau_AcoOff[i] = new TH1D(Form("FA_tau_AcoOff_%i", i), "", 120, 60, 120);
		Eta_tau_AcoOff[i] = new TH1D(Form("Eta_tau_AcoOff_%i", i), "", 120, 60, 120);
	}

	cout << "entires is " << s->t1->GetEntries() << endl;

	for (int i = 0; i < s->t1->GetEntries(); i++)
	{
		double percentage = 100.0 * (i - 0) / (s->t1->GetEntries() - 0);

		if ((i - 0) % 1000 == 0)
		{

			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = now - start;

			// Estimate remaining time based on entries processed within the specified range
			double remaining_time = elapsed.count() * (s->t1->GetEntries() - i) / (i - 0 + 1);

			// Convert remaining time to a human-readable format (hours, minutes, seconds)
			int hours = static_cast<int>(remaining_time) / 3600;
			int minutes = (static_cast<int>(remaining_time) % 3600) / 60;
			int seconds = static_cast<int>(remaining_time) % 60;

			// Output progress and estimated time remaining
			std::cout << "Progress: " << percentage << "% completed, "
					  << "Estimated time remaining: "
					  << hours << "h " << minutes << "m " << seconds << "s\r"
					  << std::flush;
		}

		// W part
		s->t1->GetEntry(i);

		// event selection cut
		// hfCoincFilter2Th4 = 1,
		// primaryVertexFilter = 2,
		// clusterCompatibilityFilter = 3,

		if (!(s->evtSel[1]))
			continue;
		if (!(s->evtSel[2]))
			continue;
		if (!(s->evtSel[3]))
			continue;
		if (abs(s->bestvtxZ) > 15)
			continue;

		int hiBin = 0;
		hiBin = s->getcentrality(3);

		float eventweight = 1.0;
		eventweight = vzRW->reweightFactor(s->bestvtxZ) * vzRW->reweightFactorCent(hiBin) * s->Ncoll[hiBin] * (s->weightLHE_gen->at(1080) / 10000.0);
		event_weight->Fill(0.5, eventweight);

		if (!(s->trigHLT[6]))
			continue;

		bool isTau = 0;
		for (unsigned int j = 0; j < s->candSize_gen; j++)
		{
			if (TMath::Abs(s->DecayID_gen[j]) == 23)
			{
				double ptWeight = spectrumRW->getReweightFactorMuon(s->pT_gen[j]);
				eventweight *= ptWeight;
				isTau = false;
				break;
			}
			if (TMath::Abs(s->DecayID_gen[j]) == 15)
			{
				isTau = 1;
				break;
			}
		}

		// This is for reco level
		for (unsigned int j = 0; j < s->candSize; j++)
		{
			if (s->mass[j] < s->masslowlimit || s->mass[j] > s->masshighlimit)
				continue;
			if (s->pTD1[j] < s->ptlowlimit)
				continue;
			if (s->pTD2[j] < s->ptlowlimit)
				continue;
			if (abs(s->EtaD1[j]) > 2.4)
				continue;
			if (abs(s->EtaD2[j]) > 2.4)
				continue;
			if (s->VtxProb[j] < 0.001)
				continue;
			if (abs(s->y[j]) > 2.4)
				continue;
			if (!(s->tightMuon1Cut(j, "POG") && s->tightMuon2Cut(j, "POG")))
				continue;

			bool isDaughter1Trigger = s->trigMuon1->at(6).at(j); // 6 is HLT_HIL3Mu12;
			bool isDaughter2Trigger = s->trigMuon2->at(6).at(j);
			if (!(isDaughter1Trigger || isDaughter2Trigger))
				continue;

			bool isOppositeSign = s->chargeD1[j] != s->chargeD2[j];

			int centbinpositioncounter[11] = {};
			s->CentBinSearching(centbinpositioncounter, hiBin);

			float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(s->PhiD1[j] - s->PhiD2[j]))) / TMath::Pi();
			bool passesAco[3] = {1, 1, 1};
			if (s->pT[j] < 1.25 && acoplanarity < 0.001)
				passesAco[0] = false;

			if (isOppositeSign)
			{
				for (int k = 0; k < s->centarraysize; k++)
				{
					if (centbinpositioncounter[k] != 0)
					{
						if (!isTau)
						{
							if (passesAco[0])
							{
								double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
								FA_nominal[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);

								if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
								{
									double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
									Eta_nominal[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);
								}
							}

							double efficiency_acooff = s->getEfficiency(e_without_aco[k], s->y[j], s->pT[j]);
							FA_AcoOff[k]->Fill(s->mass[j], 1.0 / efficiency_acooff * eventweight);

							if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
							{
								double efficiency_acooff = s->getEfficiency(e_without_aco[k], s->y[j], s->pT[j]);
								Eta_AcoOff[k]->Fill(s->mass[j], 1.0 / efficiency_acooff * eventweight);
							}
						}
						if (isTau)
						{
							if (passesAco[0])
							{
								double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
								FA_tau_nominal[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);

								if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
								{
									double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
									Eta_tau_nominal[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);
								}
							}

							double efficiency_acooff = s->getEfficiency(e_without_aco[k], s->y[j], s->pT[j]);
							FA_tau_AcoOff[k]->Fill(s->mass[j], 1.0 / efficiency_acooff * eventweight);

							if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
							{
								double efficiency_acooff = s->getEfficiency(e_without_aco[k], s->y[j], s->pT[j]);
								Eta_tau_AcoOff[k]->Fill(s->mass[j], 1.0 / efficiency_acooff * eventweight);
							}
						}
					}
				}
			}
		}
	}

	TFile *histogram_file;
	if (opt == 1)
		histogram_file = new TFile("./rootfile/mc_signal.root", "UPDATE");
	if (opt == 2)
		histogram_file = new TFile("./rootfile/mc_w.root", "UPDATE");
	if (opt == 3)
		histogram_file = new TFile("./rootfile/mc_tt.root", "UPDATE");

	histogram_file->cd(); // Skipped Oct 8
	event_weight->Write("", 2);

	for (int i = 0; i < s->centarraysize; i++)
	{
		FA_nominal[i]->Write("", 2);
		Eta_nominal[i]->Write("", 2);
		FA_AcoOff[i]->Write("", 2);
		Eta_AcoOff[i]->Write("", 2);

		FA_tau_nominal[i]->Write("", 2);
		Eta_tau_nominal[i]->Write("", 2);
		FA_tau_AcoOff[i]->Write("", 2);
		Eta_tau_AcoOff[i]->Write("", 2);
	}

	histogram_file->Close();
	s->f1->Close();
}
