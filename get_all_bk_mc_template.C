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

void get_all_bk_mc_template(int type = 4)
{
	// 1 = PbPb
	// 2 = pp bk
	// 3 = pp nobk
	// 4 = pp mass range

	// 5 = PbPb zoomin nominal
	// 6 = PbPb zoomin tnpD
	// 7 = PbPb zoomin tnpU
	// 8 = PbPb zoomin acooff
	// 9 = PbPb zoomin nominal_no_bk_sub
	// 10 = PbPb zoomin nominal_binning
	// 11 = PbPb zoomin nominal_range

	// 12 = pp zoomin nominal
	// 13 = pp zoomin tnpD
	// 14 = pp zoomin tnpU
	// 15 = pp zoomin acooff
	// 16 = pp zoomin nominal_no_bk_sub
	// 17 = pp zoomin nominal_binning
	// 18 = pp zoomin nominal_range

	MC_18 *s = new MC_18();

	auto start = std::chrono::high_resolution_clock::now();
	TH1::SetDefaultSumw2();
	// opt == 1 means signal MC;
	// opt == 2 means W;
	// opt == 3 means tt;
	// dont have to worry about starlight

	const int nbins_mass_shift = 42;
	const int nbins_smear = 42;
	const int nbins_cent = 11;

	TFile *ratio_file_raw[nbins_cent];
	TFile *ratio_file_eta[nbins_cent];

	TString postname = "";

	for (int i = 0; i < nbins_cent; i++)
	{
		if (!((i < 4) || (i == 10)))
			continue;
		if (type == 1)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/raw_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/eta_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}
		else if (type == 2)
		{
			ratio_file_raw[i] = new TFile("../bwgaus/pp_raw_bk.root", "READ");
			ratio_file_eta[i] = new TFile("../bwgaus/pp_eta_bk.root", "READ");
		}
		else if (type == 3)
		{
			ratio_file_raw[i] = new TFile("../bwgaus/pp_raw_nobk.root", "READ");
			ratio_file_eta[i] = new TFile("../bwgaus/pp_eta_nobk.root", "READ");
		}
		else if (type == 4)
		{
			ratio_file_raw[i] = new TFile("../bwgaus/pp_raw_mass_range.root", "READ");
			ratio_file_eta[i] = new TFile("../bwgaus/pp_eta_mass_range.root", "READ");
		}

		// Zoom in PbPb

		else if (type == 5)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin/PbPb_FA_nominal_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin/PbPb_eta_nominal_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}
		else if (type == 6)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin/PbPb_FA_tnpU_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin/PbPb_eta_tnpU_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}
		else if (type == 7)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin/PbPb_FA_tnpD_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin/PbPb_eta_tnpD_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}
		else if (type == 8)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin/PbPb_FA_acooff_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin/PbPb_eta_acooff_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}
		else if (type == 9)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin/PbPb_FA_nominal_no_bk_sub_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin/PbPb_eta_nominal_no_bk_sub_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}
		else if (type == 10)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin/PbPb_FA_nominal_binning_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin/PbPb_eta_nominal_binning_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}
		else if (type == 11)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin/PbPb_FA_nominal_range_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin/PbPb_eta_nominal_range_%i_%i.root", s->cenlowlimit[i], s->cenhighlimit[i]), "READ");
		}

		// Zoom in pp

		else if (type == 12)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_FA_nominal_%i.root", 22), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_eta_nominal_%i.root", 22), "READ");
		}
		else if (type == 13)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_FA_tnpU_%i.root", 22), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_eta_tnpU_%i.root", 22), "READ");
		}
		else if (type == 14)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_FA_tnpD_%i.root", 22), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_eta_tnpD_%i.root", 22), "READ");
		}
		else if (type == 15)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_FA_acooff_%i.root", 22), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_eta_acooff_%i.root", 22), "READ");
		}
		else if (type == 16)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_FA_nominal_no_bk_sub_%i.root", 22), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_eta_nominal_no_bk_sub_%i.root", 22), "READ");
		}
		else if (type == 17)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_FA_nominal_binning_%i.root", 22), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_eta_nominal_binning_%i.root", 22), "READ");
		}
		else if (type == 18)
		{
			ratio_file_raw[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_FA_nominal_range_%i.root", 22), "READ");
			ratio_file_eta[i] = new TFile(Form("../bwgaus/zoomin_pp/pp_eta_nominal_range_%i.root", 22), "READ");
		}
	}

	TString ratio_name = "Compare_shift_%i_smear_%i";
	TF1 *t_ratio_raw[nbins_mass_shift][nbins_smear][nbins_cent];
	TF1 *t_ratio_eta[nbins_mass_shift][nbins_smear][nbins_cent];

	TH1D *template_FA_nominal[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *template_FA_acooff[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *template_FA_mass_range[nbins_mass_shift][nbins_smear][nbins_cent];

	TH1D *template_Eta_nominal[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *template_Eta_acooff[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *template_Eta_mass_range[nbins_mass_shift][nbins_smear][nbins_cent];

	for (int i = 0; i < nbins_mass_shift; ++i)
	{
		for (int j = 0; j < nbins_smear; j++)
		{
			for (int k = 0; k < nbins_cent; k++)
			{
				if (!((k < 4) || (k == 10)))
					continue;
				t_ratio_raw[i][j][k] = (TF1 *)ratio_file_raw[k]->Get(Form(ratio_name, i, j));
				t_ratio_eta[i][j][k] = (TF1 *)ratio_file_eta[k]->Get(Form(ratio_name, i, j));

				template_FA_nominal[i][j][k] = new TH1D(Form("template_FA_nominal_%i_%i_%i", i, j, k), "", 120, 60, 120);
				template_FA_acooff[i][j][k] = new TH1D(Form("template_FA_acooff_%i_%i_%i", i, j, k), "", 120, 60, 120);
				template_Eta_nominal[i][j][k] = new TH1D(Form("template_Eta_nominal_%i_%i_%i", i, j, k), "", 120, 60, 120);
				template_Eta_acooff[i][j][k] = new TH1D(Form("template_Eta_acooff_%i_%i_%i", i, j, k), "", 120, 60, 120);

				template_FA_mass_range[i][j][k] = new TH1D(Form("template_FA_mass_range_%i_%i_%i", i, j, k), "", 80, 70, 110);
				template_Eta_mass_range[i][j][k] = new TH1D(Form("template_Eta_mass_range_%i_%i_%i", i, j, k), "", 80, 70, 110);
			}
		}
	}

	gSystem->Load("./header/libDict.so");

	MCReweight *vzRW = new MCReweight("WeightsAndEfficiencies_Z2mummu/vzReweight.root", "WeightsAndEfficiencies_Z2mummu/centralityFlatteningWeight.root");
	PtReweightSpectrum *spectrumRW = new PtReweightSpectrum("WeightsAndEfficiencies_Z2mummu/ptSpectrumReweighting.root");
	MuonTnP *tnp = new MuonTnP();

	s->SetupRootfile(1, 1);
	s->SetupBranches(0);

	TEfficiency *e[11];
	TEfficiency *e_acooff[11];

	TFile *eff_f1 = new TFile("./rootfile/mc_eff.root", "READ");
	for (int i = 0; i < 11; i++)
	{
		e[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
		e_acooff[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_noAco_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
	}

	// Histogram setup
	TH1D *event_weight = new TH1D("event_weight", "event_weight", 1, 0, 1);

	cout << "entires is " << s->t1->GetEntries() << endl;

	Long64_t totalEntries = s->t1->GetEntries();

	for (int i = 0; i < totalEntries; i++)
	{
		double percentage = 100.0 * (i - 0) / (totalEntries - 0);

		if ((i - 0) % 1000 == 0)
		{

			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = now - start;

			// Estimate remaining time based on entries processed within the specified range
			double remaining_time = elapsed.count() * (totalEntries - i) / (i - 0 + 1);

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
		// cout << eventweight << endl;
		event_weight->Fill(0.5, eventweight);

		if (!(s->trigHLT[6]))
			continue;

		bool isTau = 0;
		for (unsigned int j = 0; j < s->candSize_gen; j++)
		{
			if (TMath::Abs(s->DecayID_gen[j]) == 23)
			{
				double ptWeight = spectrumRW->getReweightFactorMuon(s->pT_gen[j]);
				// cout << "ptWeight is " << ptWeight << endl;
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
			// cout << "hiBin is " << hiBin << endl;

			float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(s->PhiD1[j] - s->PhiD2[j]))) / TMath::Pi();
			bool passesAco[3] = {1, 1, 1};
			if (s->pT[j] < 1.25 && acoplanarity < 0.001)
				passesAco[0] = false;

			if (isOppositeSign)
			{
				bool isgenmatching = false;
				double matchedgenindex = -99;
				for (int gen_index = 0; gen_index < s->candSize_gen; gen_index++)
				{
					if (s->RecIdx_gen[gen_index] == j)
					{
						isgenmatching = true;
						matchedgenindex = gen_index;
						break;
					}
				}

				if (!isgenmatching)
					continue;

				for (int k = 0; k < s->centarraysize; k++)
				{
					if (!((k < 4) || (k == 10)))
						continue;
					if (centbinpositioncounter[k] != 0)
					{
						if (!isTau)
						{
							if (passesAco[0])
							{
								double gen_mass = s->Calc_Z_gen(matchedgenindex);
								double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);

								for (int massshift = 0; massshift < nbins_mass_shift; massshift++)
								{
									for (int smear = 0; smear < nbins_smear; smear++)
									{
										double ratio_weight = t_ratio_raw[massshift][smear][k]->Eval(gen_mass);
										template_FA_nominal[massshift][smear][k]->Fill(s->mass[j], 1.0 / efficiency * eventweight * ratio_weight);
										template_FA_mass_range[massshift][smear][k]->Fill(s->mass[j], 1.0 / efficiency * eventweight * ratio_weight);
									}
								}

								if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
								{
									for (int massshift = 0; massshift < nbins_mass_shift; massshift++)
									{
										for (int smear = 0; smear < nbins_smear; smear++)
										{
											double ratio_weight = t_ratio_eta[massshift][smear][k]->Eval(gen_mass);
											template_Eta_nominal[massshift][smear][k]->Fill(s->mass[j], 1.0 / efficiency * eventweight * ratio_weight);
											template_Eta_mass_range[massshift][smear][k]->Fill(s->mass[j], 1.0 / efficiency * eventweight * ratio_weight);
										}
									}
								}
							}

							double gen_mass = s->Calc_Z_gen(matchedgenindex);
							double efficiency_acooff = s->getEfficiency(e_acooff[k], s->y[j], s->pT[j]);

							for (int massshift = 0; massshift < nbins_mass_shift; massshift++)
							{
								for (int smear = 0; smear < nbins_smear; smear++)
								{
									double ratio_weight = t_ratio_raw[massshift][smear][k]->Eval(gen_mass);
									template_FA_acooff[massshift][smear][k]->Fill(s->mass[j], 1.0 / efficiency_acooff * eventweight * ratio_weight);
								}
							}

							if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
							{
								for (int massshift = 0; massshift < nbins_mass_shift; massshift++)
								{
									for (int smear = 0; smear < nbins_smear; smear++)
									{
										double ratio_weight = t_ratio_eta[massshift][smear][k]->Eval(gen_mass);
										template_Eta_acooff[massshift][smear][k]->Fill(s->mass[j], 1.0 / efficiency_acooff * eventweight * ratio_weight);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	TFile *histogram_file;
	if (type == 1)
		histogram_file = new TFile("./rootfile/template_PbPb.root", "UPDATE");
	if (type == 2)
		histogram_file = new TFile("./rootfile/template_pp_bk.root", "UPDATE");
	if (type == 3)
		histogram_file = new TFile("./rootfile/template_pp_nobk.root", "UPDATE");
	if (type == 4)
		histogram_file = new TFile("./rootfile/template_pp_mass_range.root", "UPDATE");

	if (type == 5)
		histogram_file = new TFile("./rootfile/template_PbPb_zoomin_nominal.root", "UPDATE");
	if (type == 6)
		histogram_file = new TFile("./rootfile/template_PbPb_zoomin_tnpU.root", "UPDATE");
	if (type == 7)
		histogram_file = new TFile("./rootfile/template_PbPb_zoomin_tnpD.root", "UPDATE");
	if (type == 8)
		histogram_file = new TFile("./rootfile/template_PbPb_zoomin_acooff.root", "UPDATE");
	if (type == 9)
		histogram_file = new TFile("./rootfile/template_PbPb_zoomin_nominal_no_bksub.root", "UPDATE");
	if (type == 10)
		histogram_file = new TFile("./rootfile/template_PbPb_zoomin_nominal_binning.root", "UPDATE");
	if (type == 11)
		histogram_file = new TFile("./rootfile/template_PbPb_zoomin_nominal_range.root", "UPDATE");

	if (type == 12)
		histogram_file = new TFile("./rootfile/template_pp_zoomin_nominal.root", "UPDATE");
	if (type == 13)
		histogram_file = new TFile("./rootfile/template_pp_zoomin_tnpU.root", "UPDATE");
	if (type == 14)
		histogram_file = new TFile("./rootfile/template_pp_zoomin_tnpD.root", "UPDATE");
	if (type == 15)
		histogram_file = new TFile("./rootfile/template_pp_zoomin_acooff.root", "UPDATE");
	if (type == 16)
		histogram_file = new TFile("./rootfile/template_pp_zoomin_nominal_no_bksub.root", "UPDATE");
	if (type == 17)
		histogram_file = new TFile("./rootfile/template_pp_zoomin_nominal_binning.root", "UPDATE");
	if (type == 18)
		histogram_file = new TFile("./rootfile/template_pp_zoomin_nominal_range.root", "UPDATE");

	histogram_file->cd();

	for (int i = 0; i < nbins_cent; i++)
	{
		if (!((i < 4) || (i == 10)))
			continue;

		for (int j = 0; j < nbins_mass_shift; j++)
		{
			for (int k = 0; k < nbins_smear; k++)
			{
				template_FA_nominal[j][k][i]->Write("", 2);
				template_FA_acooff[j][k][i]->Write("", 2);
				template_Eta_nominal[j][k][i]->Write("", 2);
				template_Eta_acooff[j][k][i]->Write("", 2);
				template_FA_mass_range[j][k][i]->Write("", 2);
				template_Eta_mass_range[j][k][i]->Write("", 2);
			}
		}
	}
}
