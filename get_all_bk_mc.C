#include "MC_18.h"
void get_all_bk_mc(){
	TH1::SetDefaultSumw2();
	MC_18 *w_bk_obj = new MC_18();
	MC_18 *tt_bk_obj = new MC_18();
	MC_18 *mc_signal_obj = new MC_18(); // tau is here

	tt_bk_obj->SetupRootfile(3,1);
	tt_bk_obj->SetupBranches(0);

	w_bk_obj->SetupRootfile(2,1);
	w_bk_obj->SetupBranches(0);

	mc_signal_obj->SetupRootfile(1,1);
	mc_signal_obj->SetupBranches(0);

	//Histogram setup
	//W area
	TH1D *W_event_weight = new TH1D("W_event_weight","W_event_weight",1,0,1);
	TH1D *mass_array_W[11];
	for (int i = 0; i < tt_bk_obj->centarraysize; i++){
		mass_array_W[i] = new TH1D(Form("mass_array_W_%i",i),Form("mc_W_bk_%i",w_bk_obj->cenhighlimit[i]),120,60,120);
	}
	//TT area
	TH1D *tt_event_weight = new TH1D("tt_event_weight","tt_event_weight",1,0,1);
	TH1D *mass_array_tt[11];
	for (int i = 0; i < tt_bk_obj->centarraysize; i++){
		mass_array_tt[i] = new TH1D(Form("mass_array_tt_%i",i),Form("mc_tt_bk_%i",tt_bk_obj->cenhighlimit[i]),120,60,120);
	}
	//DY area
	TH1D *signal_event_weight = new TH1D("signal_event_weight","signal_event_weight",1,0,1);
	TH1D *mass_array_tau[11];
	TH1D *mass_array_signal[11];
	TH1D *mc_reco_A[11];
	for (int i = 0; i < tt_bk_obj->centarraysize; i++){
		mass_array_tau[i] = new TH1D(Form("mass_array_tau_%i",i),Form("mc_tau_bk_%i",mc_signal_obj->cenhighlimit[i]),120,60,120);
		mass_array_signal[i] = new TH1D(Form("mass_array_signal_%i",i),Form("mc_signal_bk_%i",mc_signal_obj->cenhighlimit[i]),120,60,120);
		mc_reco_A[i] = new TH1D(Form("mc_reco_A_%i",i),Form("mc_reco_A_%i",mc_signal_obj->cenhighlimit[i]),1000,0,1);
	}

	cout << "TT entires is " << tt_bk_obj->t1->GetEntries() << endl;
	cout << "W entires is " << w_bk_obj->t1->GetEntries() << endl;
	cout << "DY entires is " << mc_signal_obj->t1->GetEntries() << endl;
	cout << "W has the most entries, use this to start looping" << endl;

	for (int i = 0; i < w_bk_obj->t1->GetEntries(); i++){
		if (i%100000 == 0) cout << "Event " << i << endl;

		// W part
		w_bk_obj->t1->GetEntry(i);
		W_event_weight->Fill(0.5,w_bk_obj->weightLHE_gen->at(1080)/10000.0);
		bool isTau_W = 0;
		for (int j=0; j<w_bk_obj->candSize_gen; j++){
			if (w_bk_obj->DecayID_gen[j] == 15) {
				isTau_W = 1;
				break;
			}
		}
		if (!isTau_W){
			for (int j=0;j<w_bk_obj->candSize; j++){
				if(w_bk_obj->mass[j] < w_bk_obj->masslowlimit || w_bk_obj->mass[j] > w_bk_obj->masshighlimit) continue;
				if(w_bk_obj->pTD1[j] < w_bk_obj->ptlowlimit) continue;
				if(w_bk_obj->pTD2[j] < w_bk_obj->ptlowlimit) continue;
				if(abs(w_bk_obj->EtaD1[j]) > 2.4)continue;
				if(abs(w_bk_obj->EtaD2[j]) > 2.4)continue;
				if(w_bk_obj->VtxProb[j] < 0.001)continue;
			
				int centbinpositioncounter[11] = {};
				w_bk_obj->CentBinSearching(centbinpositioncounter);
				for (int k=0; k<w_bk_obj->centarraysize; k++){
					if (centbinpositioncounter[k] != 0 ){
						mass_array_W[k]->Fill(w_bk_obj->mass[j],w_bk_obj->weightLHE_gen->at(1080)/10000.0);
					}
				}
			}
		}

		// TT part
		if (i < tt_bk_obj->t1->GetEntries()){
			tt_bk_obj->t1->GetEntry(i);
			tt_event_weight->Fill(0.5,tt_bk_obj->weightLHE_gen->at(1080)/10000.0);

			bool isTau_tt = 0;
			for (int j=0; j<tt_bk_obj->candSize_gen; j++){
				if (tt_bk_obj->DecayID_gen[j] == 15) {
					isTau_tt = 1;
					break;
				}
			}
			if (!isTau_tt){
				for (int j=0;j<tt_bk_obj->candSize; j++){
					if(tt_bk_obj->mass[j] < tt_bk_obj->masslowlimit || tt_bk_obj->mass[j] > tt_bk_obj->masshighlimit) continue;
					if(tt_bk_obj->pTD1[j] < tt_bk_obj->ptlowlimit) continue;
					if(tt_bk_obj->pTD2[j] < tt_bk_obj->ptlowlimit) continue;
					if(abs(tt_bk_obj->EtaD1[j]) > 2.4)continue;
					if(abs(tt_bk_obj->EtaD2[j]) > 2.4)continue;
					if(tt_bk_obj->VtxProb[j] < 0.001)continue;

					int centbinpositioncounter[11] = {};
					tt_bk_obj->CentBinSearching(centbinpositioncounter);
					for (int k=0; k<tt_bk_obj->centarraysize; k++){
						if (centbinpositioncounter[k] != 0 ){
							mass_array_tt[k]->Fill(tt_bk_obj->mass[j],tt_bk_obj->weightLHE_gen->at(1080)/10000.0);
						}
					}
				}
			}
		}
		// DY part
		if (i < mc_signal_obj->t1->GetEntries()){
			mc_signal_obj->t1->GetEntry(i);
			signal_event_weight->Fill(0.5,mc_signal_obj->weightLHE_gen->at(1080)/10000.0);
			bool isTau_signal = 0;
			for (int j=0; j<mc_signal_obj->candSize_gen; j++){
				if (mc_signal_obj->DecayID_gen[j] == 15) {
					isTau_signal = 1;
					break;
				}
			}
			if (isTau_signal){
				for (int j=0;j<mc_signal_obj->candSize; j++){
					if(mc_signal_obj->mass[j] < mc_signal_obj->masslowlimit || mc_signal_obj->mass[j] > mc_signal_obj->masshighlimit) continue;
					if(mc_signal_obj->pTD1[j] < mc_signal_obj->ptlowlimit) continue;
					if(mc_signal_obj->pTD2[j] < mc_signal_obj->ptlowlimit) continue;
					if(abs(mc_signal_obj->EtaD1[j]) > 2.4)continue;
					if(abs(mc_signal_obj->EtaD2[j]) > 2.4)continue;
					if(mc_signal_obj->VtxProb[j] < 0.001)continue;

					int centbinpositioncounter[11] = {};
					mc_signal_obj->CentBinSearching(centbinpositioncounter);
					for (int k=0; k<mc_signal_obj->centarraysize; k++){
						if (centbinpositioncounter[k] != 0){
							mass_array_tau[k]->Fill(mc_signal_obj->mass[j],mc_signal_obj->weightLHE_gen->at(1080)/10000.0);
						}
					}
				}
			}
			else{
				for (int j=0;j<mc_signal_obj->candSize; j++){
					if(mc_signal_obj->mass[j] < mc_signal_obj->masslowlimit || mc_signal_obj->mass[j] > mc_signal_obj->masshighlimit) continue;
					if(mc_signal_obj->pTD1[j] < mc_signal_obj->ptlowlimit) continue;
					if(mc_signal_obj->pTD2[j] < mc_signal_obj->ptlowlimit) continue;
					if(abs(mc_signal_obj->EtaD1[j]) > 2.4)continue;
					if(abs(mc_signal_obj->EtaD2[j]) > 2.4)continue;
					if(mc_signal_obj->VtxProb[j] < 0.001)continue;

					int centbinpositioncounter[11] = {};
					mc_signal_obj->CentBinSearching(centbinpositioncounter);
					for (int k=0; k<mc_signal_obj->centarraysize; k++){
						if (centbinpositioncounter[k] != 0){
							mass_array_signal[k]->Fill(mc_signal_obj->mass[j],mc_signal_obj->weightLHE_gen->at(1080)/10000.0);

							if ((abs(mc_signal_obj->PhiD1[j]) + abs(mc_signal_obj->PhiD2[j])) > 3.1415926){
								mc_reco_A[k]->Fill(1-((2*3.1415926 - abs(mc_signal_obj->PhiD1[j]) - abs(mc_signal_obj->PhiD2[j]))/ 3.1415926));
							}
							else{
								mc_reco_A[k]->Fill(1-((abs(mc_signal_obj->PhiD1[j]) + abs(mc_signal_obj->PhiD2[j]))/ 3.1415926));
							}
						}
					}
				}
			}
		}
	}

	TFile *histogram_file = new TFile("./rootfile/mc_bk.root","UPDATE");
	histogram_file->cd();
	W_event_weight->Write("",2);
	tt_event_weight->Write("",2);
	signal_event_weight->Write("",2);
	for (int i = 0; i < mc_signal_obj->centarraysize; i++){
		mass_array_W[i]->Write("",2);
		mass_array_tt[i]->Write("",2);
		mass_array_tau[i]->Write("",2);
		mass_array_signal[i]->Write("",2);
		mc_reco_A[i]->Write("",2);
	}
	histogram_file->Close();
	tt_bk_obj->f1->Close();
	w_bk_obj->f1->Close();
	mc_signal_obj->f1->Close();
}