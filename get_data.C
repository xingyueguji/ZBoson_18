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
	gSystem->Load("./libDict.so");

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

	TFile *eff_f1 = new TFile("./rootfile/mc_eff.root", "READ");
	for (int i = 0; i < 11; i++)
	{
		e[i] = (TEfficiency *)eff_f1->Get(Form("eff_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]));
	}

	// TH1D *data_A[11];
	TH1D *mass_array_data[11];
	TH1D *mass_array_data_withy[11];
	TH1D *mass_array_data_witheta[11];
	TH1D *mass_array_data_with_eff[11];
	TH1D *mass_array_data_withy_witheff[11];
	TH1D *mass_array_data_witheta_witheff[11];

	TH1D *mass_array_data_same_sign[11];
	TH1D *mass_array_data_same_sign_withy[11];
	TH1D *mass_array_data_same_sign_witheta[11];
	TH1D *mass_array_data_same_sign_with_eff[11];
	TH1D *mass_array_data_same_sign_withy_witheff[11];
	TH1D *mass_array_data_same_sign_witheta_witheff[11];

	TH1D *y_without_cut;
	TH1D *y_with_cut;

	TH1D *h_cent = new TH1D("h_cent", "h_cent", 100, 0, 200);

	TH1D *h_Z_pt = new TH1D("h_Z_pt", "h_Z_pt", 300, 0, 300);
	TH1D *h_Z_pt_eta = new TH1D("h_Z_pt_eta", "h_Z_pt_eta", 300, 0, 300);

	y_with_cut = new TH1D("y_with_cut", "Z_rapidity_with_cut", 100, -3, 3);
	y_without_cut = new TH1D("y_without_cut", "Z_rapidity_without_cut", 100, -3, 3);

	// TH1D *h_pesudorapidity = new TH1D("pesudorapidity","",100,-10,10);
	// TH1F *h_rapidity = new TH1F("rapidity","",16,-2.4,2.4);

	for (int i = 0; i < data->centarraysize; i++)
	{
		mass_array_data[i] = new TH1D(Form("mass_array_data_%i", i), Form("mass_data_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_withy[i] = new TH1D(Form("mass_array_data_withy_%i", i), Form("mass_data_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_with_eff[i] = new TH1D(Form("mass_array_data_with_eff_%i", i), Form("mass_data_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_withy_witheff[i] = new TH1D(Form("mass_array_data_withy_witheff_%i", i), Form("mass_data_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_witheta[i] = new TH1D(Form("mass_array_data_witheta_%i", i), Form("mass_data_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_witheta_witheff[i] = new TH1D(Form("mass_array_data_witheta_witheff_%i", i), Form("mass_data_%i_%i", data->cenlowlimit[i], data->cenhighlimit[i]), 120, 60, 120);

		mass_array_data_same_sign[i] = new TH1D(Form("mass_array_data_same_sign_%i", i), Form("mass_data_same_sign_%i_%i", data_same_sign->cenlowlimit[i], data_same_sign->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_same_sign_withy[i] = new TH1D(Form("mass_array_data_same_sign_withy_%i", i), Form("mass_data_same_sign_%i_%i", data_same_sign->cenlowlimit[i], data_same_sign->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_same_sign_with_eff[i] = new TH1D(Form("mass_array_data_same_sign_with_eff_%i", i), Form("mass_data_same_sign_%i_%i", data_same_sign->cenlowlimit[i], data_same_sign->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_same_sign_withy_witheff[i] = new TH1D(Form("mass_array_data_same_sign_withy_witheff_%i", i), Form("mass_data_same_sign_%i_%i", data_same_sign->cenlowlimit[i], data_same_sign->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_same_sign_witheta[i] = new TH1D(Form("mass_array_data_same_sign_witheta_%i", i), Form("mass_data_same_sign_%i_%i", data_same_sign->cenlowlimit[i], data_same_sign->cenhighlimit[i]), 120, 60, 120);
		mass_array_data_same_sign_witheta_witheff[i] = new TH1D(Form("mass_array_data_same_sign_witheta_witheff_%i", i), Form("mass_data_same_sign_%i_%i", data_same_sign->cenlowlimit[i], data_same_sign->cenhighlimit[i]), 120, 60, 120);
	}

	// Place for the new unbinned dataset:
	RooRealVar *roomass[11];
	RooDataSet *roodata[11];
	RooDataSet *roodata_y[11];
	RooDataSet *roodata_eff[11];
	RooDataSet *roodata_y_eff[11];
	RooDataSet *roodata_eta[11];
	RooDataSet *roodata_eta_eff[11];
	RooDataSet *roodata_raw_uniform[11];

	for (int i = 0; i < data->centarraysize; i++)
	{
		roomass[i] = new RooRealVar("roomass", "roomass", 60, 120);
		roomass[i]->setBinning(RooBinning(10000, 60, 120));
		roodata[i] = new RooDataSet(Form("roodata_%i", i), Form("roodata_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		roodata_y[i] = new RooDataSet(Form("roodata_y_%i", i), Form("roodata_y_%i", i), RooArgSet(*roomass[i]));
		roodata_eff[i] = new RooDataSet(Form("roodata_eff_%i", i), Form("roodata_eff_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		roodata_y_eff[i] = new RooDataSet(Form("roodata_y_eff_%i", i), Form("roodata_y_eff_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		roodata_eta[i] = new RooDataSet(Form("roodata_eta_%i", i), Form("roodata_eta_%i", i), RooArgSet(*roomass[i]));
		roodata_eta_eff[i] = new RooDataSet(Form("roodata_eta_eff_%i", i), Form("roodata_eta_eff_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
		roodata_raw_uniform[i] = new RooDataSet(Form("roodata_raw_uniform_%i", i), Form("roodata_raw_uniform_%i", i), RooArgSet(*roomass[i]), RooFit::WeightVar("eventWeight"));
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
		hiBin = data->getcentrality(0); // FIXME: doZDC? hiBinVar?

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

			float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(data->PhiD1[j] - data->PhiD2[j]))) / TMath::Pi();
			bool passesAco[3] = {1, 1, 1};
			if (data->pT[j] < 1.25 && acoplanarity < 0.001)
				passesAco[0] = false;

			h_cent->Fill(hiBin);

			// This is the place for eta cut check:

			y_without_cut->Fill(data->y[j]);
			if ((abs(data->EtaD1[j]) < 1) && (abs(data->EtaD2[j]) < 1))
			{
				if (abs(data->y[j]) > 1)
					cout << "Warning Z has y > 1 with eta cut < 1 on both" << endl;
				y_with_cut->Fill(data->y[j]);
				h_Z_pt_eta->Fill(data->pT[j]);
			}
			h_Z_pt->Fill(data->pT[j]);

			for (int k = 0; k < data->centarraysize; k++)
			{
				if (centbinpositioncounter[k] != 0)
				{
					// with y cut
					if (passesAco[0])
					{
						if (abs(data->y[j]) < 1)
						{
							double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
							mass_array_data_withy[k]->Fill(data->mass[j]);
							mass_array_data_withy_witheff[k]->Fill(data->mass[j], 1.0 / efficiency);
							roomass[k]->setVal(data->mass[j]);
							roodata_y[k]->add(RooArgSet(*roomass[k]));
							roodata_y_eff[k]->add(RooArgSet(*roomass[k]), 1.0 / efficiency);
						}
						// with eta cut < 1 on both D1 and D2 muons
						if ((abs(data->EtaD1[j]) < 1) && (abs(data->EtaD2[j]) < 1))
						{
							double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
							mass_array_data_witheta[k]->Fill(data->mass[j]);
							mass_array_data_witheta_witheff[k]->Fill(data->mass[j], 1.0 / efficiency);
							roomass[k]->setVal(data->mass[j]);
							roodata_eta[k]->add(RooArgSet(*roomass[k]));
							roodata_eta_eff[k]->add(RooArgSet(*roomass[k]), 1.0 / efficiency);
						}
						// without y cut
						double efficiency = data->getEfficiency(e[k], data->y[j], data->pT[j]);
						mass_array_data[k]->Fill(data->mass[j]);
						mass_array_data_with_eff[k]->Fill(data->mass[j], 1.0 / efficiency);

						// This is the raw one: test validity of unbin fit
						roomass[k]->setVal(data->mass[j]);
						roodata[k]->add(RooArgSet(*roomass[k]));
						roodata_eff[k]->add(RooArgSet(*roomass[k]), 1.0 / efficiency);
						roodata_raw_uniform[k]->add(RooArgSet(*roomass[k]), 1);
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
		hiBin = data_same_sign->getcentrality(0); // FIXME: doZDC? hiBinVar?

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
					// with y cut
					if (passesAco[0])
					{
						if (abs(data_same_sign->y[j]) < 1)
						{
							double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
							mass_array_data_same_sign_withy[k]->Fill(data_same_sign->mass[j]);
							mass_array_data_same_sign_withy_witheff[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
						}
						// with eta cut
						if ((abs(data_same_sign->EtaD1[j]) < 1) && (abs(data_same_sign->EtaD2[j]) < 1))
						{
							double efficiency = data_same_sign->getEfficiency(e[k], data_same_sign->y[j], data_same_sign->pT[j]);
							mass_array_data_same_sign_witheta[k]->Fill(data_same_sign->mass[j]);
							mass_array_data_same_sign_witheta_witheff[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
						}
						// without y cut
						double efficiency = data_same_sign->getEfficiency(e[k], data->y[j], data->pT[j]);
						mass_array_data_same_sign[k]->Fill(data_same_sign->mass[j]);
						mass_array_data_same_sign_with_eff[k]->Fill(data_same_sign->mass[j], 1.0 / efficiency);
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
		mass_array_data_same_sign_withy[i]->Write("", 2);
		mass_array_data_same_sign_withy_witheff[i]->Write("", 2);
		mass_array_data_same_sign[i]->Write("", 2);
		mass_array_data_same_sign_with_eff[i]->Write("", 2);
		mass_array_data_same_sign_witheta[i]->Write("", 2);
		mass_array_data_same_sign_witheta_witheff[i]->Write("", 2);

		mass_array_data_withy[i]->Write("", 2);
		mass_array_data_withy_witheff[i]->Write("", 2);
		mass_array_data[i]->Write("", 2);
		mass_array_data_with_eff[i]->Write("", 2);
		mass_array_data_witheta[i]->Write("", 2);
		mass_array_data_witheta_witheff[i]->Write("", 2);

		roodata[i]->Write("", 2);
		roodata_y[i]->Write("", 2);
		roodata_y_eff[i]->Write("", 2);
		roodata_eta[i]->Write("", 2);
		roodata_eta_eff[i]->Write("", 2);
		roodata_eff[i]->Write("", 2);
		roodata_raw_uniform[i]->Write("", 2);
	}

	TCanvas *c1 = new TCanvas("", "", 1200, 600);
	c1->Divide(2, 2);
	// c1->SetLogy(1);
	c1->cd(1);
	h_Z_pt_eta->Draw("hist");
	c1->cd(2);
	h_Z_pt->Draw("hist");
	c1->cd(3);
	y_with_cut->Draw("hist");
	c1->cd(4);
	y_without_cut->Draw("hist");

	TCanvas *c2 = new TCanvas("", "", 1200, 600);

	c2->cd();
	h_cent->Draw();

	c1->SaveAs("./etacheck/data_pt.pdf");
	c2->SaveAs("./etacheck/cent.pdf");

	histogram_file->Close();
	data->f1->Close();
}