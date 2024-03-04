#include "MC_18.h"

void MC_BK(Int_t opt = 2,Int_t iteration = 0){
	TH1::SetDefaultSumw2();
	MC_18 *obj1 = new MC_18();
	obj1->SetupRootfile(opt);
	obj1->SetupBranches();

	//variable declaration
	Int_t itlow = 32;
	Int_t ithigh = 34;

	//Histogram
	TH1D *mc_reco_Zmass = new TH1D("mc_reco_Zmass","mc_reco_Zmass",120,60,120);
	TH1D *h_p_1 = new TH1D("h_p_1","Pull",120,60,120);
	TH1D *mc_gen_Zmass = new TH1D("mc_gen_Zmass","mc_gen_Zmass",120,61.1876,121.1876);
	TH1D *h_p_2 = new TH1D("h_p_2","Pull",120,61.1876,121.1876);

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

	if (iteration == 0) itlow = 0, ithigh = nentries;
	else{
		itlow = (iteration - 1) * 100000, ithigh = iteration * 100000;
		if (iteration == 189) ithigh = nentries;
		//else itlow = (iteration - 1) * 150000 + 1, ithigh = iteration * 150000;
		cout << "itlow is " << itlow << " ithigh is " << ithigh << endl;
	}

	for (int i= itlow; i<ithigh; i++){
		if (i%100000 == 0) cout << "Event " << i << endl;
		//cout << "Event " << i << endl;
		obj1->t1->GetEntry(i);

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
			//cout <<"PIDD1 is " << obj1->PIDD1[j] << " PIDD2 is " << obj1->PIDD2[j] << endl;
			// Why PIDD1 or PIDD2 == -77 ???
			//if (abs(obj1->PIDD1[j]) != 13 ||obj1->PIDD2[j] != 13) cout << "PID for reco D1 is " << obj1->PIDD1[j] << " PID for reco D2 is " << obj1->PIDD2[j] << endl;
			mc_reco_Zmass->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);

		}
	}

	//Draw plot
	TFile *file = new TFile("./mc_bk_plot/W_histogram.root","UPDATE");
	file->cd();
	mc_reco_Zmass->Write(Form("mc_W_bk_%i",iteration),2);
	file->Close();

}