#include "MC_18.h"
void MC_BK_tt(){
	TH1::SetDefaultSumw2();
	MC_18 *obj1 = new MC_18();
	obj1->SetupRootfile(3,1);
	obj1->SetupBranches(2);

	//variable declaration
	
	//Histogram
	TH1D *mc_reco_Zmass = new TH1D("mc_reco_Zmass","mc_reco_Zmass",120,60,120);
	TH1D *mc_reco_Zmass_0_24 = new TH1D("mc_reco_Zmass_0_24","mc_reco_Zmass_0_24",120,60,120);
	TH1D *mc_reco_Zmass_24_100 = new TH1D("mc_reco_Zmass_24_100","mc_reco_Zmass_24_100",120,60,120);

	TH1D *nEvents = new TH1D("nEvents","nEvents",1,0,1);

	//Canvas

	//TPad

	//Event loop
	int nentries = obj1->t1->GetEntries();
	cout << "We have " << nentries << " events in total" << endl;

	for (int i= 0; i<nentries; i++){
		if (i%100000 == 0) cout << "Event " << i << endl;
		//cout << "Event " << i << endl;
		obj1->t1->GetEntry(i);

		nEvents->Fill(0.5,obj1->weightLHE_gen->at(1080)/10000.0);

		//Cut setup		
		bool isTau = 0;
		for (int k=0; k<obj1->candSize_gen; k++){
			//cout << "PID gen is " << obj1->PID_gen[k] << "PID D1 is " << obj1->PIDD1_gen[k] << " PID D2 is " << obj1->PIDD2_gen[k]<< endl;
			//if (abs(obj1->PID_gen[k]) == 24) cout << "This is W event" << endl;
			// Why PID_gen is 221 mostly?
			if (obj1->DecayID_gen[k] == 15) {
				//cout << "Skip Tau Event !" << endl;
				isTau = 1;
				break;
			}
		}
		if (isTau) continue;

		for (int j=0;j<obj1->candSize; j++){
			if(obj1->mass[j] < obj1->masslowlimit || obj1->mass[j] > obj1->masshighlimit) continue;
			if(obj1->pTD1[j] < obj1->ptlowlimit) continue;
			if(obj1->pTD2[j] < obj1->ptlowlimit) continue;
			if(abs(obj1->EtaD1[j]) > 2.4)continue;
			if(abs(obj1->EtaD2[j]) > 2.4)continue;
			if(obj1->VtxProb[j] < 0.001)continue;
			if (obj1->centrality/2 >= 0 && obj1->centrality/2 <24) mc_reco_Zmass_0_24->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
			if (obj1->centrality/2 >= 24 && obj1->centrality/2 <100)mc_reco_Zmass_24_100->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
			//cout <<"PIDD1 is " << obj1->PIDD1[j] << " PIDD2 is " << obj1->PIDD2[j] << endl;
			// Why PIDD1 or PIDD2 == -77 ???
			//if (abs(obj1->PIDD1[j]) != 13 ||obj1->PIDD2[j] != 13) cout << "PID for reco D1 is " << obj1->PIDD1[j] << " PID for reco D2 is " << obj1->PIDD2[j] << endl;
			mc_reco_Zmass->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);

		}
	}

	//Draw plot
	TFile *file = new TFile("./mc_bk_plot/W_histogram.root","UPDATE");
	file->cd();
	mc_reco_Zmass->Write("mc_bk_tt",2);
	mc_reco_Zmass_0_24->Write("mc_bk_tt_24",2);
	mc_reco_Zmass_24_100->Write("mc_bk_tt_100",2);
	nEvents->Write("tt_weight_factor",2);
	file->Close();
	obj1->f1->Close();
}