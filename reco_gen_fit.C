#include "MC_18.h"
void reco_gen_fit(bool userapiditycut = 0){
	MC_18 *mc_signal_obj = new MC_18(); // tau is here
	mc_signal_obj->SetupRootfile(1,1);
	mc_signal_obj->SetupBranches(0);

	TH1D *mc_reco_hist = new TH1D("mc_reco_hist","",120,60,120);
	TH1D *mc_gen_hist = new TH1D("mc_gen_hist","",120,60,120);

	for (int i = 0; i < mc_signal_obj->t1->GetEntries();i++){
			mc_signal_obj->t1->GetEntry(i);
			if (i % 100000 == 0) cout << i << endl;
			bool isTau_signal = 0;
			mc_signal_obj->Calc_Z_gen(mc_gen_hist);
			for (int j=0; j<mc_signal_obj->candSize_gen; j++){
				if (mc_signal_obj->DecayID_gen[j] == 15) {
					isTau_signal = 1;
					break;
				}
			}

			if (!isTau_signal){
				for (int j=0;j<mc_signal_obj->candSize; j++){
					if(mc_signal_obj->mass[j] < mc_signal_obj->masslowlimit || mc_signal_obj->mass[j] > mc_signal_obj->masshighlimit) continue;
					if(mc_signal_obj->pTD1[j] < mc_signal_obj->ptlowlimit) continue;
					if(mc_signal_obj->pTD2[j] < mc_signal_obj->ptlowlimit) continue;
					if(abs(mc_signal_obj->EtaD1[j]) > 2.4)continue;
					if(abs(mc_signal_obj->EtaD2[j]) > 2.4)continue;
					if(mc_signal_obj->VtxProb[j] < 0.001)continue;
					//new Z rapidity cut:
					if (userapiditycut){
						mc_signal_obj->getrapidity(j);
						if(abs(mc_signal_obj->rapidity) > 1)continue;
					}
					mc_reco_hist->Fill(mc_signal_obj->mass[j],mc_signal_obj->weightLHE_gen->at(1080)/10000.0);
				}
			}
		}
	TFile *f12 = new TFile("./justplot.root","UPDATE");
	f12->cd();
	mc_reco_hist->Write("",2);
	mc_gen_hist->Write("",2);

	TCanvas *c1 = new TCanvas("","",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	mc_reco_hist->Draw("HIST");
	mc_gen_hist->Draw("HIST");
	c1->SaveAs("justplot.pdf");
}