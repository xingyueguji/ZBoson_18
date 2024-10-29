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

void get_all_bk_mc(int jobID = 0, int numSegments = 100, int opt = 1, int numberofsamples = 1000, double shiftlowbin = -0.5, double shifthighbin = 0.2, double smearlowbin = 0.0, double smearhighbin = 0.015)
{

	auto start = std::chrono::high_resolution_clock::now();
	TH1::SetDefaultSumw2();
	// opt == 1 means signal MC;
	// opt == 2 means W;
	// opt == 3 means tt;
	// dont have to worry about starlight

	// Create matrix of shifted/smeared signal mc mass distribution.
	const int nbins_mass_shift = 21;
	const int nbins_smear = 21;
	const int nbins_cent = 11;
	cout << "This is mass shift from " << shiftlowbin << " to " << shifthighbin << " Smearing from " << smearlowbin << " to " << smearhighbin << " With " << numberofsamples << " Samples " << " Binning : " << nbins_mass_shift << " " << nbins_smear << endl;
	TH1D *modifiedmass[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *modifiedmass_raw[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *modifiedmass_raw_without_eff[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *modifiedmass_eta_without_eff[nbins_mass_shift][nbins_smear][nbins_cent];

	TH1D *modifiedmass_raw_without_eff_new[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *modifiedmass_eta_without_eff_new[nbins_mass_shift][nbins_smear][nbins_cent];

	for (int i = 0; i < nbins_mass_shift; i++)
	{
		for (int j = 0; j < nbins_smear; j++)
		{
			for (int k = 0; k < nbins_cent; k++)
			{
				modifiedmass[i][j][k] = new TH1D(Form("modifiedmass_%i_%i_%i", i, j, k), "", 120, 60, 120);
				modifiedmass_raw[i][j][k] = new TH1D(Form("modifiedmass_raw_%i_%i_%i", i, j, k), "", 120, 60, 120);
				modifiedmass_raw_without_eff[i][j][k] = new TH1D(Form("modifiedmass_raw_without_eff_%i_%i_%i", i, j, k), "", 120, 60, 120);
				modifiedmass_eta_without_eff[i][j][k] = new TH1D(Form("modifiedmass_eta_without_eff_%i_%i_%i", i, j, k), "", 120, 60, 120);
				modifiedmass_raw_without_eff_new[i][j][k] = new TH1D(Form("modifiedmass_raw_without_eff_%i_%i_%i_new", i, j, k), "", 120, 60, 120);
				modifiedmass_eta_without_eff_new[i][j][k] = new TH1D(Form("modifiedmass_eta_without_eff_%i_%i_%i_new", i, j, k), "", 120, 60, 120);
			}
		}
	}

	double shift_bin[nbins_mass_shift];
	double smear_bin[nbins_smear];

	double shiftamount = (abs(shiftlowbin) + abs(shifthighbin)) / (nbins_mass_shift - 1); // 21 ticks from -0.5 to 0.5, 20 divisions
	double smearedamount = (smearhighbin) / (nbins_smear - 1);							  // 21 ticks from 0 to 0.2, 20 divisions.

	for (int i = 0; i < nbins_mass_shift; i++)
	{
		shift_bin[i] = shiftlowbin + i * shiftamount;
	}
	for (int j = 0; j < nbins_smear; j++)
	{
		smear_bin[j] = 0 + j * smearedamount;
	}

	TRandom3 randGen;

	double mean = 1;

	gSystem->Load("./libDict.so");

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

	TFile *eff_f1 = new TFile("./rootfile/mc_eff.root", "READ");
	for (int i = 0; i < 11; i++)
	{
		e[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
	}

	// Histogram setup
	TH1D *event_weight = new TH1D("event_weight", "event_weight", 1, 0, 1);
	TH1D *mass_array[11];
	TH1D *mass_array_withy[11];
	TH1D *mass_array_with_eff[11];
	TH1D *mass_array_withy_witheff[11];
	TH1D *mass_array_witheta[11];
	TH1D *mass_array_witheta_witheff[11];

	TH1D *mass_array_tau[11];
	TH1D *mass_array_tau_withy[11];
	TH1D *mass_array_tau_with_eff[11];
	TH1D *mass_array_tau_withy_witheff[11];
	TH1D *mass_array_tau_witheta[11];
	TH1D *mass_array_tau_witheta_witheff[11];

	TH1D *y_without_cut;
	TH1D *y_with_cut;

	RooRealVar *roomass[11];

	RooDataSet *rooreco[11];
	RooDataSet *rooreco_y[11];
	RooDataSet *rooreco_eff[11];
	RooDataSet *rooreco_y_eff[11];
	RooDataSet *rooreco_eta[11];
	RooDataSet *rooreco_eta_eff[11];

	RooDataSet *roogen[11];
	RooDataSet *roogen_y[11];
	RooDataSet *roogen_eta[11];

	for (int i = 0; i < s->centarraysize; i++)
	{
		roomass[i] = new RooRealVar("roomass", "roomass", 60, 120);
		roomass[i]->setBinning(RooBinning(10000, 60, 120));

		rooreco[i] = new RooDataSet(Form("rooreco_%i", i), Form("rooreco_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		rooreco_y[i] = new RooDataSet(Form("rooreco_y_%i", i), Form("rooreco_y_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		rooreco_eff[i] = new RooDataSet(Form("rooreco_eff_%i", i), Form("rooreco_eff_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		rooreco_y_eff[i] = new RooDataSet(Form("rooreco_y_eff_%i", i), Form("rooreco_y_eff_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		rooreco_eta[i] = new RooDataSet(Form("rooreco_eta_%i", i), Form("rooreco_eta_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		rooreco_eta_eff[i] = new RooDataSet(Form("rooreco_eta_eff_%i", i), Form("rooreco_eta_eff_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));

		roogen[i] = new RooDataSet(Form("roogen_%i", i), Form("roogen_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		roogen_y[i] = new RooDataSet(Form("roogen_y_%i", i), Form("roogen_y_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		roogen_eta[i] = new RooDataSet(Form("roogen_eta_%i", i), Form("roogen_eta_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
	}

	y_with_cut = new TH1D("y_with_cut", "Z_rapidity_with_cut", 100, -3, 3);
	y_without_cut = new TH1D("y_without_cut", "Z_rapidity_without_cut", 100, -3, 3);

	for (int i = 0; i < s->centarraysize; i++)
	{
		mass_array[i] = new TH1D(Form("mass_array_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_withy[i] = new TH1D(Form("mass_array_withy_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_with_eff[i] = new TH1D(Form("mass_array_with_eff_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_withy_witheff[i] = new TH1D(Form("mass_array_withy_witheff_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_witheta[i] = new TH1D(Form("mass_array_witheta_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_witheta_witheff[i] = new TH1D(Form("mass_array_witheta_witheff_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);

		mass_array_tau[i] = new TH1D(Form("mass_array_tau_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_tau_withy[i] = new TH1D(Form("mass_array_tau_withy_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_tau_with_eff[i] = new TH1D(Form("mass_array_tau_with_eff_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_tau_withy_witheff[i] = new TH1D(Form("mass_array_tau_withy_witheff_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_tau_witheta[i] = new TH1D(Form("mass_array_tau_witheta_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
		mass_array_tau_witheta_witheff[i] = new TH1D(Form("mass_array_tau_witheta_witheff_%i", i), Form("mc_bk_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]), 120, 60, 120);
	}

	cout << "entires is " << s->t1->GetEntries() << endl;

	Long64_t totalEntries = s->t1->GetEntries();
	Long64_t entriesPerSegment = totalEntries / numSegments;
	Long64_t startEntry = jobID * entriesPerSegment;
	Long64_t endEntry = (jobID == numSegments - 1) ? totalEntries : (jobID + 1) * entriesPerSegment;
	std::cout << "Processing segment " << jobID << ": Entries " << startEntry << " to " << endEntry - 1 << std::endl;

	for (int i = startEntry; i < endEntry; i++)
	{
		double percentage = 100.0 * i / s->t1->GetEntries();
		if (i % 1000 == 0)
		{
			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = now - start;

			// Estimate remaining time
			double remaining_time = elapsed.count() * (s->t1->GetEntries() - i) / i;

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

		// This is for gen level, skipped oct 8 just for faster processing
		/*for (unsigned int j = 0; j < s->candSize_gen; j++)
		{
			if (s->PID_gen[j] != 23)
				continue;
			if (s->DecayID_gen[j] != 23)
				continue;
			if (s->pTD1_gen[j] < s->ptlowlimit)
				continue;
			if (s->pTD2_gen[j] < s->ptlowlimit)
				continue;
			if (abs(s->EtaD1_gen[j]) > 2.4)
				continue;
			if (abs(s->EtaD2_gen[j]) > 2.4)
				continue;
			if (abs(s->y_gen[j]) > 2.4)
				continue;

			bool isOppositeSign = s->chargeD1_gen[j] != s->chargeD2_gen[j];
			int centbinpositioncounter[11] = {};
			s->CentBinSearching(centbinpositioncounter, hiBin);

			if (isOppositeSign)
			{
				for (int k = 0; k < s->centarraysize; k++)
				{
					if (centbinpositioncounter[k] != 0)
					{
						Float_t massofZ_gen = s->Calc_Z_gen(j);
						if (massofZ_gen >= 120 || massofZ_gen < 60)
							continue;
						if (abs(s->y_gen[j]) < 1)
						{
							roomass[k]->setVal(massofZ_gen);
							roogen_y[k]->add(RooArgSet(*roomass[k]), eventweight);
						}
						if ((abs(s->EtaD1_gen[j]) < 1) && (abs(s->EtaD2_gen[j]) < 1))
						{
							roomass[k]->setVal(massofZ_gen);
							roogen_eta[k]->add(RooArgSet(*roomass[k]), eventweight);
						}
						// raw one
						roomass[k]->setVal(massofZ_gen);
						roogen[k]->add(RooArgSet(*roomass[k]), eventweight);
					}
				}
			}
		}*/

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

			if (isOppositeSign)
			{
				// Here's for eta distribution check skipped Oct 8
				/*if (opt == 1)
				{
					y_without_cut->Fill(s->y[j]);
					if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
					{
						if (abs(s->y[j]) > 1)
							cout << "Warning Z has y > 1 with eta cut < 1 on both" << endl;
						y_with_cut->Fill(s->y[j]);
					}
				}*/

				for (int k = 0; k < s->centarraysize; k++)
				{
					if (centbinpositioncounter[k] != 0)
					{
						if (!isTau)
						{
							/*double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
							mass_array[k]->Fill(s->mass[j], eventweight);
							mass_array_with_eff[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);

							roomass[k]->setVal(s->mass[j]);
							rooreco[k]->add(RooArgSet(*roomass[k]), eventweight);
							rooreco_eff[k]->add(RooArgSet(*roomass[k]), 1.0 / efficiency * eventweight);*/

							for (int l = 0; l < nbins_mass_shift; l++)
							{
								for (int m = 0; m < nbins_smear; m++)
								{
									for (int nsamples = 0; nsamples < numberofsamples; ++nsamples)
									{
										double randomsmear = randGen.Gaus(mean, smear_bin[m]);
										double updatedmass = (s->mass[j]) * randomsmear + shift_bin[l];

										double randomsmear_1 = randGen.Gaus(0, smear_bin[m]) * 91.1876;
										double updatedmass_1 = (s->mass[j]) + randomsmear_1 + shift_bin[l];
										// cout << "updatedmass is " << updatedmass << endl;
										//  modifiedmass_raw[l][m][k]->Fill(updatedmass, 1.0 / efficiency * eventweight);
										modifiedmass_raw_without_eff[l][m][k]->Fill(updatedmass, eventweight * (1.0 / numberofsamples));
										modifiedmass_raw_without_eff_new[l][m][k]->Fill(updatedmass_1, eventweight * (1.0 / numberofsamples));
									}
								}
							}

							if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
							{
								// Skipped Oct 8
								/*double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
								mass_array_witheta[k]->Fill(s->mass[j], eventweight);
								mass_array_witheta_witheff[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);

								roomass[k]->setVal(s->mass[j]);
								rooreco_eta[k]->add(RooArgSet(*roomass[k]), eventweight);
								rooreco_eta_eff[k]->add(RooArgSet(*roomass[k]), 1.0 / efficiency * eventweight);*/

								// Here's the shift and smearing part:
								for (int l = 0; l < nbins_mass_shift; l++)
								{
									for (int m = 0; m < nbins_smear; m++)
									{
										for (int nsamples = 0; nsamples < numberofsamples; ++nsamples)
										{
											double randomsmear = randGen.Gaus(mean, smear_bin[m]);
											double updatedmass = (s->mass[j]) * randomsmear + shift_bin[l];
											double randomsmear_1 = randGen.Gaus(0, smear_bin[m]) * 91.1876;
											double updatedmass_1 = (s->mass[j]) + randomsmear_1 + shift_bin[l];
											// modifiedmass[l][m][k]->Fill(updatedmass, 1.0 / efficiency * eventweight);
											modifiedmass_eta_without_eff[l][m][k]->Fill(updatedmass, eventweight * (1.0 / numberofsamples));
											modifiedmass_eta_without_eff_new[l][m][k]->Fill(updatedmass_1, eventweight * (1.0 / numberofsamples));
										}
									}
								}
							}
							// Skipped Oct 8
							/*if (abs(s->y[j]) < 1)
							{
								mass_array_withy[k]->Fill(s->mass[j], eventweight);
								mass_array_withy_witheff[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);

								roomass[k]->setVal(s->mass[j]);
								rooreco_y[k]->add(RooArgSet(*roomass[k]), eventweight);
								rooreco_y_eff[k]->add(RooArgSet(*roomass[k]), 1.0 / efficiency * eventweight);
							}*/
						}
						/*if (isTau) // Skipped Oct 8
						{
							double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
							mass_array_tau[k]->Fill(s->mass[j], eventweight);
							mass_array_tau_with_eff[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);

							if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1))
							{
								double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
								mass_array_tau_witheta[k]->Fill(s->mass[j], eventweight);
								mass_array_tau_witheta_witheff[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);
							}

							if (abs(s->y[j]) < 1)
							{
								mass_array_tau_withy[k]->Fill(s->mass[j], eventweight);
								mass_array_tau_withy_witheff[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);
							}
						}*/
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

	/*histogram_file->cd(); // Skipped Oct 8
	event_weight->Write("", 2);

	for (int i = 0; i < s->centarraysize; i++)
	{
		mass_array[i]->Write("", 2);
		mass_array_withy[i]->Write("", 2);
		mass_array_with_eff[i]->Write("", 2);
		mass_array_withy_witheff[i]->Write("", 2);
		mass_array_witheta[i]->Write("", 2);
		mass_array_witheta_witheff[i]->Write("", 2);

		mass_array_tau[i]->Write("", 2);
		mass_array_tau_withy[i]->Write("", 2);
		mass_array_tau_with_eff[i]->Write("", 2);
		mass_array_tau_withy_witheff[i]->Write("", 2);
		mass_array_tau_witheta[i]->Write("", 2);
		mass_array_tau_witheta_witheff[i]->Write("", 2);

		rooreco[i]->Write("", 2);
		rooreco_y[i]->Write("", 2);
		rooreco_eff[i]->Write("", 2);
		rooreco_y_eff[i]->Write("", 2);
		rooreco_eta[i]->Write("", 2);
		rooreco_eta_eff[i]->Write("", 2);

		roogen[i]->Write("", 2);
		roogen_y[i]->Write("", 2);
		roogen_eta[i]->Write("", 2);
	}*/

	TFile *modified_histogram_file = new TFile(Form("./rootfile/jobID_%i_shift_%.1f_%.1f_smear_%.1f_%.3f_modified_signal_%i_%i_%i.root", jobID, shiftlowbin ,shifthighbin, smearlowbin, smearhighbin, nbins_mass_shift, nbins_smear, numberofsamples), "UPDATE");
	modified_histogram_file->cd();

	for (int i = 0; i < nbins_mass_shift; i++)
	{
		for (int j = 0; j < nbins_smear; j++)
		{
			for (int k = 0; k < nbins_cent; k++)
			{
				// modifiedmass[i][j][k]->Write("", 2);
				// modifiedmass_raw[i][j][k]->Write("", 2);
				modifiedmass_raw_without_eff[i][j][k]->Write("", 2);
				modifiedmass_eta_without_eff[i][j][k]->Write("", 2);

				modifiedmass_raw_without_eff_new[i][j][k]->Write("", 2);
				modifiedmass_eta_without_eff_new[i][j][k]->Write("", 2);
			}
		}
	}

	/*if (opt == 1) // Skipped Oct 8
	{
		TCanvas *c1 = new TCanvas("", "", 1200, 600);
		c1->Divide(2, 1);
		// c1->SetLogy(1);
		c1->cd(1);
		y_without_cut->Draw("hist");
		c1->cd(2);
		y_with_cut->Draw("hist");
		c1->SaveAs("./etacheck/mc_signal.pdf");
	}*/

	histogram_file->Close();
	modified_histogram_file->Close();
	s->f1->Close();

	cout << "All Finished F*** the internet" << endl;
}
