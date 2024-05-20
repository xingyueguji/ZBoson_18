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

void get_all_bk_mc(int opt = 3){
	TH1::SetDefaultSumw2();
	//opt == 1 means signal MC;
	//opt == 2 means W;
	//opt == 3 means tt;
	//dont have to worry about starlight
	gSystem->Load("./libDict.so");

	MC_18 *s = new MC_18();

	MCReweight *vzRW = new MCReweight("WeightsAndEfficiencies_Z2mummu/vzReweight.root","WeightsAndEfficiencies_Z2mummu/centralityFlatteningWeight.root");
	PtReweightSpectrum *spectrumRW = new PtReweightSpectrum("WeightsAndEfficiencies_Z2mummu/ptSpectrumReweighting.root");
	MuonTnP *tnp = new MuonTnP();
	

	if (opt == 1) {
		s->SetupRootfile(1,1);
		s->SetupBranches(0);
	}
	if (opt == 2){
		s->SetupRootfile(2,1);
		s->SetupBranches(0);
	}
	if (opt == 3){
		s->SetupRootfile(3,1);
		s->SetupBranches(0);
	}

	TEfficiency* e[11];

	TFile *eff_f1 = new TFile("./rootfile/mc_eff.root","READ");
	for (int i=0; i<11; i++){
		e[i] = (TEfficiency*)eff_f1->Get(Form("eff_noSF_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]));
	}

	//Histogram setup
	TH1D *event_weight = new TH1D("event_weight","event_weight",1,0,1);
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

	y_with_cut = new TH1D("y_with_cut","Z_rapidity_with_cut",100,-3,3);
	y_without_cut = new TH1D("y_without_cut","Z_rapidity_without_cut",100,-3,3);

	for (int i = 0; i < s->centarraysize; i++){
		mass_array[i] = new TH1D(Form("mass_array_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_withy[i] = new TH1D(Form("mass_array_withy_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_with_eff[i] = new TH1D(Form("mass_array_with_eff_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_withy_witheff[i] = new TH1D(Form("mass_array_withy_witheff_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_witheta[i]= new TH1D(Form("mass_array_witheta_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_witheta_witheff[i]= new TH1D(Form("mass_array_witheta_witheff_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);

		mass_array_tau[i] = new TH1D(Form("mass_array_tau_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_tau_withy[i] = new TH1D(Form("mass_array_tau_withy_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_tau_with_eff[i] = new TH1D(Form("mass_array_tau_with_eff_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_tau_withy_witheff[i] = new TH1D(Form("mass_array_tau_withy_witheff_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_tau_witheta[i] = new TH1D(Form("mass_array_tau_witheta_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
		mass_array_tau_witheta_witheff[i] = new TH1D(Form("mass_array_tau_witheta_witheff_%i",i),Form("mc_bk_%i_%i",s->cenlowlimit[i],s->cenhighlimit[i]),120,60,120);
	}

	cout << "entires is " << s->t1->GetEntries() << endl;

	for (int i = 0; i < s->t1->GetEntries(); i++){
		if (i%100000 == 0) cout << "Event " << i << endl;

		// W part
		s->t1->GetEntry(i);

		//event selection cut
		//hfCoincFilter2Th4 = 1,
        //primaryVertexFilter = 2,
        //clusterCompatibilityFilter = 3,
		if (!(s->evtSel[1])) continue;
		if (!(s->evtSel[2])) continue;
		if (!(s->evtSel[3])) continue;
		if (abs(s->bestvtxZ) > 15) continue;

		int hiBin = 0;
		hiBin = s->getcentrality(3);
		
		float eventweight = 1.0;
		eventweight = vzRW->reweightFactor(s->bestvtxZ) * vzRW->reweightFactorCent(hiBin) * s->Ncoll[hiBin] * (s->weightLHE_gen->at(1080)/10000.0);
		event_weight->Fill(0.5,eventweight);

		if (!(s->trigHLT[6])) continue;

		bool isTau = 0;
		for (unsigned int j=0; j<s->candSize_gen; j++){
			if (TMath::Abs(s->DecayID_gen[j]) == 23) {
				double ptWeight = spectrumRW->getReweightFactorMuon(s->pT_gen[j]);
				eventweight *= ptWeight;
				isTau = false;
				break;
			}
			if (TMath::Abs(s->DecayID_gen[j]) == 15) {
				isTau = 1;
				break;
			}
		}
		for (unsigned int j=0;j<s->candSize; j++){
			if(s->mass[j] < s->masslowlimit || s->mass[j] > s->masshighlimit) continue;
			if(s->pTD1[j] < s->ptlowlimit) continue;
			if(s->pTD2[j] < s->ptlowlimit) continue;
			if(abs(s->EtaD1[j]) > 2.4)continue;
			if(abs(s->EtaD2[j]) > 2.4)continue;
			if(s->VtxProb[j] < 0.001)continue;
			if (abs(s->y[j]) > 2.4) continue;
			if (!(s->tightMuon1Cut(j,"POG") && s->tightMuon2Cut(j,"POG"))) continue;

			bool isDaughter1Trigger = s->trigMuon1->at(6).at(j); // 6 is HLT_HIL3Mu12;
			bool isDaughter2Trigger = s->trigMuon2->at(6).at(j);
			if (!(isDaughter1Trigger || isDaughter2Trigger)) continue;

			bool isOppositeSign =  s->chargeD1[j] != s->chargeD2[j];

			int centbinpositioncounter[11] = {};
			s->CentBinSearching(centbinpositioncounter,hiBin);
			
			if (isOppositeSign){
				//Here's for eta distribution check
				if (opt == 1){
					y_without_cut->Fill(s->y[j]);
					if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j]) < 1)) {
						if (abs(s->y[j]) > 1) cout << "Warning Z has y > 1 with eta cut < 1 on both" << endl;
						y_with_cut->Fill(s->y[j]);
					}	
				}

				for (int k = 0; k < s->centarraysize; k++){
					if (centbinpositioncounter[k] != 0){
						if (!isTau){
							double efficiency = s->getEfficiency(e[k],s->y[j],s->pT[j]);
							mass_array[k]->Fill(s->mass[j],eventweight);
							mass_array_with_eff[k]->Fill(s->mass[j],1.0/efficiency * eventweight);

							if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j])< 1)){
								double efficiency = s->getEfficiency(e[k],s->y[j],s->pT[j]);
								mass_array_witheta[k]->Fill(s->mass[j],eventweight);
								mass_array_witheta_witheff[k]->Fill(s->mass[j],1.0/efficiency * eventweight);
							}

							if (abs(s->y[j]) < 1){
								mass_array_withy[k]->Fill(s->mass[j],eventweight);
								mass_array_withy_witheff[k]->Fill(s->mass[j],1.0/efficiency * eventweight);
							}
						}
						if (isTau){
							double efficiency = s->getEfficiency(e[k],s->y[j],s->pT[j]);
							mass_array_tau[k]->Fill(s->mass[j],eventweight);
							mass_array_tau_with_eff[k]->Fill(s->mass[j],1.0/efficiency * eventweight);

							if ((abs(s->EtaD1[j]) < 1) && (abs(s->EtaD2[j])< 1)){
								double efficiency = s->getEfficiency(e[k],s->y[j],s->pT[j]);
								mass_array_tau_witheta[k]->Fill(s->mass[j],eventweight);
								mass_array_tau_witheta_witheff[k]->Fill(s->mass[j],1.0/efficiency * eventweight);
							}

							if (abs(s->y[j]) < 1){
								mass_array_tau_withy[k]->Fill(s->mass[j],eventweight);
								mass_array_tau_withy_witheff[k]->Fill(s->mass[j],1.0/efficiency * eventweight);
							}
						}
					}				
				}
			}
		}
	}

	TFile *histogram_file;
	if (opt == 1) histogram_file = new TFile("./rootfile/mc_signal.root","UPDATE");
	if (opt == 2) histogram_file = new TFile("./rootfile/mc_w.root","UPDATE");
	if (opt == 3) histogram_file = new TFile("./rootfile/mc_tt.root","UPDATE");

	histogram_file->cd();
	event_weight->Write("",2);
	for (int i = 0; i < s->centarraysize; i++){
		mass_array[i]->Write("",2);
		mass_array_withy[i]->Write("",2);
		mass_array_with_eff[i]->Write("",2);
		mass_array_withy_witheff[i]->Write("",2);
		mass_array_witheta[i]->Write("",2);
		mass_array_witheta_witheff[i]->Write("",2);

		mass_array_tau[i]->Write("",2);
		mass_array_tau_withy[i]->Write("",2);
		mass_array_tau_with_eff[i]->Write("",2);
		mass_array_tau_withy_witheff[i]->Write("",2);
		mass_array_tau_witheta[i]->Write("",2);
		mass_array_tau_witheta_witheff[i]->Write("",2);
	}

	if (opt == 1){
		TCanvas *c1 = new TCanvas("","",600,600);
		c1->Divide(2,1);
		//c1->SetLogy(1);
		c1->cd(1);
		y_without_cut->Draw("hist");
		c1->cd(2);
		y_with_cut->Draw("hist");

		c1->SaveAs("./etacheck/mc_signal.pdf");
	}

	histogram_file->Close();
	s->f1->Close();
}