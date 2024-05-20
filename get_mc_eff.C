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

void get_mc_eff(){
	gSystem->Load("./libDict.so");

	MC_18 *mc_signal_obj = new MC_18();
	mc_signal_obj->SetupRootfile(1,1);
	mc_signal_obj->SetupBranches(0);
	MuonTnP *tnp = new MuonTnP();
	MCReweight *vzRW = new MCReweight("WeightsAndEfficiencies_Z2mummu/vzReweight.root","WeightsAndEfficiencies_Z2mummu/centralityFlatteningWeight.root");
	PtReweightSpectrum *spectrumRW = new PtReweightSpectrum("WeightsAndEfficiencies_Z2mummu/ptSpectrumReweighting.root");
	
	

	TH1::SetDefaultSumw2();

	const int centarraysize = mc_signal_obj->centarraysize;

	TH2D *recoEff_pass[centarraysize];
	TH2D *recoEff_all[centarraysize];
	TH2D *recoEff[centarraysize];

	TH2D *recoEff_pass_noSF[centarraysize];
	TH2D *recoEff_all_noSF[centarraysize];
	TH2D *recoEff_noSF[centarraysize];

	TEfficiency *eff[11];
	TEfficiency *eff_noSF[11];

	const int nrapiditybins = 16;
	const int nptbins = 9;
	int hiBin = 0;
	double ptbin[nptbins] = {0,1.0,3.0,5.0,10.0,20.0,40.0,70.0,200};

	for (unsigned int i=0 ; i<centarraysize; i++){
		recoEff_pass[i] = new TH2D(Form("recoEff_pass_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]),"",nrapiditybins,-2.4,2.4,nptbins-1,ptbin);
		recoEff_all[i] = new TH2D(Form("recoEff_all_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]),"",nrapiditybins,-2.4,2.4,nptbins-1,ptbin);
		recoEff_pass_noSF[i] = new TH2D(Form("recoEff_pass_noSF_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]),"",nrapiditybins,-2.4,2.4,nptbins-1,ptbin);
		recoEff_all_noSF[i] = new TH2D(Form("recoEff_all_noSF_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]),"",nrapiditybins,-2.4,2.4,nptbins-1,ptbin);
		//recoEff[i] = new TH2D(Form("recoEff_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]),"",nrapiditybins,-2.4,2.4,nptbins,ptbin);
	}

	for (unsigned int i=0; i < mc_signal_obj->t1->GetEntries(); i++){
		mc_signal_obj->t1->GetEntry(i);
		if (i%100000 == 0) cout << i << "/" << mc_signal_obj->t1->GetEntries() << endl;

		//event selection cut
		//hfCoincFilter2Th4 = 1,
        //primaryVertexFilter = 2,
        //clusterCompatibilityFilter = 3,
		if (!(mc_signal_obj->evtSel[1])) continue;
		if (!(mc_signal_obj->evtSel[2])) continue;
		if (!(mc_signal_obj->evtSel[3])) continue;
		if (abs(mc_signal_obj->bestvtxZ) > 15) continue;

		hiBin = mc_signal_obj->getcentrality(3);
		double eventweight = vzRW->reweightFactor(mc_signal_obj->bestvtxZ)* vzRW->reweightFactorCent(hiBin) * mc_signal_obj->Ncoll[hiBin] * mc_signal_obj->weightLHE_gen->at(1080)/10000; // Missing Cent reweight + Z reweight

		for (unsigned int j=0; j<mc_signal_obj->candSize_gen; j++){

			if (mc_signal_obj->PID_gen[j]!= 23) continue;
			if (mc_signal_obj->DecayID_gen[j]!= 23) continue;
			if (abs(mc_signal_obj->y_gen[j]) > mc_signal_obj->etalimit) continue;

			double ptWeight = spectrumRW->getReweightFactorMuon(mc_signal_obj->pT_gen[j]);
			eventweight *= ptWeight;

			if (abs(mc_signal_obj->EtaD1_gen[j]) > mc_signal_obj->etalimit) continue;
			if (abs(mc_signal_obj->EtaD2_gen[j]) > mc_signal_obj->etalimit) continue;
			if (mc_signal_obj->pTD1_gen[j] < mc_signal_obj->ptlowlimit) continue;
			if (mc_signal_obj->pTD2_gen[j] < mc_signal_obj->ptlowlimit) continue;
			int centbinpositioncounter[11] = {};
			mc_signal_obj->CentBinSearching(centbinpositioncounter,hiBin);
			for (int k=0; k<mc_signal_obj->centarraysize; k++){
				if (centbinpositioncounter[k] != 0 ){
					recoEff_all[k]->Fill(mc_signal_obj->y_gen[j],mc_signal_obj->pT_gen[j],eventweight);
					recoEff_all_noSF[k]->Fill(mc_signal_obj->y_gen[j],mc_signal_obj->pT_gen[j],eventweight);
				}
			}

			if (mc_signal_obj->RecIdx_gen[j] < 0) continue;
			int indx = mc_signal_obj->RecIdx_gen[j];

			if (mc_signal_obj->mass[indx] < mc_signal_obj->masslowlimit || mc_signal_obj->mass[indx] > mc_signal_obj->masshighlimit) continue;
			if (mc_signal_obj->pTD1[indx] < mc_signal_obj->ptlowlimit) continue;
			if (mc_signal_obj->pTD2[indx] < mc_signal_obj->ptlowlimit) continue;
			if (abs(mc_signal_obj->EtaD1[indx]) > mc_signal_obj->etalimit) continue;
			if (abs(mc_signal_obj->EtaD2[indx]) > mc_signal_obj->etalimit) continue;
			if (!(mc_signal_obj->tightMuon1Cut(indx,"POG") && mc_signal_obj->tightMuon2Cut(indx,"POG"))) continue;
			if(mc_signal_obj->VtxProb[indx] < 0.001)continue;
			// Missing Tight muon cut

			float acoplanarity =1 - TMath::Abs(TMath::ACos(TMath::Cos( mc_signal_obj->PhiD1[indx] - mc_signal_obj->PhiD2[indx] )))/TMath::Pi();
			bool passesAco[3] = {1,1,1};
			if (mc_signal_obj->pT[indx] < 1.25 && acoplanarity < 0.001) passesAco[0] = false;

			bool isDaughter1Trigger = mc_signal_obj->trigMuon1->at(6).at(indx); // 6 is HLT_HIL3Mu12;
			bool isDaughter2Trigger = mc_signal_obj->trigMuon2->at(6).at(indx);

			//cout << "trigMuon1[6][indx] is " << mc_signal_obj->trigMuon1->at(6).at(indx) << endl;

			if (!(isDaughter1Trigger || isDaughter2Trigger)) continue;

			bool isOppositeSign = mc_signal_obj->chargeD1[indx] != mc_signal_obj->chargeD2[indx];
			if (!isOppositeSign) continue;

			float scaleFactor = tnp->getZSF(mc_signal_obj->pTD1_gen[j],mc_signal_obj->EtaD1_gen[j],mc_signal_obj->pTD2_gen[j],mc_signal_obj->EtaD2_gen[j],(int)(hiBin/2),0);

			for (int k=0; k<mc_signal_obj->centarraysize; k++){
				if (centbinpositioncounter[k] != 0) {
					if (mc_signal_obj->pT_gen[j] < ptbin[8] && abs(mc_signal_obj->y_gen[j]) < 2.4){
						if (passesAco[0]){
							recoEff_pass[k]->Fill(mc_signal_obj->y_gen[j],mc_signal_obj->pT_gen[j],eventweight*scaleFactor);
							recoEff_pass_noSF[k]->Fill( mc_signal_obj->y_gen[j], mc_signal_obj->pT_gen[j],eventweight);
						}
					}
				} 
			}
		}
	}

	for (int i=0; i<centarraysize; i++){
		cout << "Now checking the consistency of histogram" << " " << i << endl;  
		mc_signal_obj->forceConsistency(recoEff_pass[i],recoEff_all[i]);
		mc_signal_obj->forceConsistency(recoEff_pass_noSF[i],recoEff_all_noSF[i]);
		recoEff[i] = (TH2D*) recoEff_pass[i]->Clone(Form("recoEff_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]));
		recoEff_noSF[i] = (TH2D*) recoEff_pass_noSF[i]->Clone(Form("recoEff_noSF_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]));
		recoEff[i]->Divide(recoEff_all[i]);
		recoEff_noSF[i]->Divide(recoEff_all_noSF[i]);

		if (TEfficiency::CheckConsistency(*(recoEff_pass[i]),*(recoEff_all[i]),"w")){
			eff[i] = new TEfficiency(*(recoEff_pass[i]),*(recoEff_all[i]));
			eff[i]->SetStatisticOption(TEfficiency::kBJeffrey);
			eff[i]->SetName(Form("eff_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]));
		}

		if (TEfficiency::CheckConsistency(*(recoEff_pass_noSF[i]),*(recoEff_all_noSF[i]),"w")){
			eff_noSF[i] = new TEfficiency(*(recoEff_pass_noSF[i]),*(recoEff_all_noSF[i]));
			eff_noSF[i]->SetStatisticOption(TEfficiency::kBJeffrey);
			eff_noSF[i]->SetName(Form("eff_noSF_%i_%i",mc_signal_obj->cenlowlimit[i],mc_signal_obj->cenhighlimit[i]));
		}

	}

	TFile *f1 = new TFile("./rootfile/mc_eff.root","UPDATE");
	f1->cd();
	for (int i=0; i<centarraysize; i++){
		recoEff[i]->Write("",2);
		recoEff_pass[i]->Write("",2);
		recoEff_all[i]->Write("",2);
		recoEff_noSF[i]->Write("",2);
		recoEff_pass_noSF[i]->Write("",2);
		recoEff_all_noSF[i]->Write("",2);

		eff[i]->Write("",2);
		eff_noSF[i]->Write("",2);
	}


}
