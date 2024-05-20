#include "MC_18.h"

void MC_plot(){
	TH1::SetDefaultSumw2();
	MC_18 *obj1 = new MC_18();
	obj1->SetupRootfile(1,1);
	obj1->SetupBranches(0);

	//variable declaration
	bool eventnotpassed = 0; // 1 means not passed, 0 means passed cut

	//Histogram
	TH1D *mc_reco_Zmass = new TH1D("mc_reco_Zmass","mc_reco_Zmass",120,60,120);
	TH1D *mc_reco_tau_bk = new TH1D("mc_reco_tau_bk","mc_reco_tau_bk",120,60,120);
	TH1D *h_p_1 = new TH1D("h_p_1","Pull",120,60,120);
	TH1D *mc_gen_Zmass = new TH1D("mc_gen_Zmass","mc_gen_Zmass",120,60,120);
	TH1D *h_p_2 = new TH1D("h_p_2","Pull",120,60,120);
	TH1D *mc_reco_Zmass_tau_0_24 = new TH1D("mc_reco_Zmass_tau_0_24","mc_reco_Zmass_tau_0_24",120,60,120);
	TH1D *mc_reco_Zmass_tau_24_100 = new TH1D("mc_reco_Zmass_tau_24_100","mc_reco_Zmass_tau_24_100",120,60,120);
	TH1D *mc_reco_Zmass_0_24 = new TH1D("mc_reco_Zmass_0_24","mc_reco_Zmass_0_24",120,60,120);
	TH1D *mc_reco_Zmass_24_100 = new TH1D("mc_reco_Zmass_24_100","mc_reco_Zmass_24_100",120,60,120);

	TH1D *mc_reco_A_plot_all = new TH1D("mc_reco_A_plot_all","mc_reco_A_plot_all",1000,0,1);
	TH1D *mc_reco_A_plot_24 = new TH1D("mc_reco_A_plot_24","mc_reco_A_plot_24",1000,0,1);
	TH1D *mc_reco_A_plot_100 = new TH1D("mc_reco_A_plot_100","mc_reco_A_plot_100",1000,0,1);

	TH1D *nEvents = new TH1D("nEvents","nEvents",1,0,1);

	//Canvas
	TCanvas *c1 = new TCanvas("","",800,600);
	TCanvas *c2 = new TCanvas("","",800,600);

	//TPad
	TPad *p_1 = new TPad("p_1","Large",0,0.3,1,1);
	TPad *p_2 = new TPad("p_2","Small",0,0,1,0.3);
	TPad *p_3 = new TPad("p_3","Large",0,0.3,1,1);
	TPad *p_4 = new TPad("p_4","Small",0,0,1,0.3);

	//Event loop
	int nentries = obj1->t1->GetEntries();
	cout << "We have " << nentries << " events in total" << endl;
	int counter = 0;

	for (int i=0; i<nentries; i++){
		if (i%100000 == 0) cout << "Event " << i << endl;
		obj1->t1->GetEntry(i);

		//if (obj1->candSize > obj1->candSize_gen) cout << "Larger" << endl;
		//obj1->Calc_Z_gen(mc_gen_Zmass); // ahead of everything so gen cut is not affect by reco cut
		if (obj1->candSize_gen == 0) cout << "NO MUON GEN" << endl;

		bool isTau = 0;
		for (int k=0; k<obj1->candSize_gen; k++){

			// Get Weight factor
			nEvents->Fill(0.5,obj1->weightLHE_gen->at(1080)/10000.0); // this is problematic !!!!
			counter++;

			if (obj1->DecayID_gen[k] == 15) {
				//cout << "Skip Tau Event !" << endl;
				isTau = 1;
				break;
			}
		}
		/*if (isTau){
			for (int j=0;j<obj1->candSize; j++){
				if(obj1->mass[j] < obj1->masslowlimit || obj1->mass[j] > obj1->masshighlimit) continue;
				if(obj1->pTD1[j] < obj1->ptlowlimit) continue;
				if(obj1->pTD2[j] < obj1->ptlowlimit) continue;
				if(abs(obj1->EtaD1[j]) > 2.4)continue;
				if(abs(obj1->EtaD2[j]) > 2.4)continue;
				if(obj1->VtxProb[j] < 0.001)continue;
				if (obj1->centrality/2 >= 0 && obj1->centrality/2 <24) mc_reco_Zmass_tau_0_24->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
			    if (obj1->centrality/2 >= 24 && obj1->centrality/2 <100)mc_reco_Zmass_tau_24_100->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
				mc_reco_tau_bk->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
			}	
		}
		else{
			for (int j=0;j<obj1->candSize; j++){
				if(obj1->mass[j] < obj1->masslowlimit || obj1->mass[j] > obj1->masshighlimit) continue;
				if(obj1->pTD1[j] < obj1->ptlowlimit) continue;
				if(obj1->pTD2[j] < obj1->ptlowlimit) continue;
				if(abs(obj1->EtaD1[j]) > 2.4)continue;
				if(abs(obj1->EtaD2[j]) > 2.4)continue;
				if(obj1->VtxProb[j] < 0.001)continue;
				if (obj1->chargeD1[j] == obj1->chargeD2[j]) cout << "SAME CHARGE EVENT WARNING!!!" << endl;
				if (obj1->centrality/2 >= 0 && obj1->centrality/2 <24){
					mc_reco_Zmass_0_24->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);

					if ((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j])) > 3.1415926){
						mc_reco_A_plot_24->Fill(1-((2*3.1415926 - abs(obj1->PhiD1[j]) - abs(obj1->PhiD2[j]))/ 3.1415926));
					}
					else{
						mc_reco_A_plot_24->Fill(1-((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j]))/ 3.1415926));
					}

				}

			    if (obj1->centrality/2 >= 24 && obj1->centrality/2 <100) {
					mc_reco_Zmass_24_100->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);

					if ((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j])) > 3.1415926){
						mc_reco_A_plot_100->Fill(1-((2*3.1415926 - abs(obj1->PhiD1[j]) - abs(obj1->PhiD2[j]))/ 3.1415926));
					}
					else{
						mc_reco_A_plot_100->Fill(1-((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j]))/ 3.1415926));
					}
				}

				mc_reco_Zmass->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
				if ((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j])) > 3.1415926){
					mc_reco_A_plot_all->Fill(1-((2*3.1415926 - abs(obj1->PhiD1[j]) - abs(obj1->PhiD2[j]))/ 3.1415926));
				}
				else{
					mc_reco_A_plot_all->Fill(1-((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j]))/ 3.1415926));
				}
			}
		}*/
	}

	cout << counter << endl;

	//Draw plot
	TFile *file = new TFile("./mc_bk_plot/W_histogram.root","UPDATE");
	file->cd();
	nEvents->Write("DY_weight_factor",2);
	/*mc_reco_Zmass->Write("ZTOUU",2);
	mc_reco_tau_bk->Write("Tau_BK",2);
	mc_gen_Zmass->Write("ZGEN",2);
	mc_reco_Zmass_0_24->Write("ZTOUU_24",2);
	mc_reco_Zmass_24_100->Write("ZTOUU_100",2);
	mc_reco_Zmass_tau_0_24->Write("Tau_BK_24",2);
	mc_reco_Zmass_tau_24_100->Write("Tau_BK_100",2);
	nEvents->Write("DY_weight_factor",2);
	mc_reco_A_plot_24->Write("MC_A_24",2);
	mc_reco_A_plot_100->Write("MC_A_100",2);
	mc_reco_A_plot_all->Write("MC_A_ALL",2);*/
	file->Close();
	obj1->f1->Close();


}