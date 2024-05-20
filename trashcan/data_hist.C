#include "MC_18.h"
void data_hist(double cendiv = 14){
	TH1::SetDefaultSumw2();
	MC_18 *obj1 = new MC_18();
	obj1->SetupRootfile(4,0);
	obj1->SetupBranches(1);

	//variable declaration
	
	//Histogram
	TH1D *mc_reco_Zmass = new TH1D("mc_reco_Zmass","mc_reco_Zmass",120,60,120);
	TH1D *mc_reco_Zmass_0_24 = new TH1D("mc_reco_Zmass_0_24","mc_reco_Zmass_0_24",120,60,120);
	TH1D *mc_reco_Zmass_24_100 = new TH1D("mc_reco_Zmass_24_100","mc_reco_Zmass_24_100",120,60,120);

	TH1D *mc_reco_A_plot_all = new TH1D("mc_reco_A_plot_all","mc_reco_A_plot_all",1000,0,1);
	TH1D *mc_reco_A_plot_24 = new TH1D("mc_reco_A_plot_24","mc_reco_A_plot_24",1000,0,1);
	TH1D *mc_reco_A_plot_100 = new TH1D("mc_reco_A_plot_100","mc_reco_A_plot_100",1000,0,1);

	//Canvas

	//TPad

	//Event loop
	int nentries = obj1->t1->GetEntries();
	cout << "We have " << nentries << " events in total" << endl;

	for (int i= 0; i<nentries; i++){
		if (i%100000 == 0) cout << "Event " << i << endl;
		//cout << "Event " << i << endl;
		obj1->t1->GetEntry(i);

		//Cut setup		

		for (int j=0;j<obj1->candSize; j++){
			if(obj1->mass[j] < obj1->masslowlimit || obj1->mass[j] > obj1->masshighlimit) continue;
			if(obj1->pTD1[j] < obj1->ptlowlimit) continue;
			if(obj1->pTD2[j] < obj1->ptlowlimit) continue;
			if(abs(obj1->EtaD1[j]) > 2.4)continue;
			if(abs(obj1->EtaD2[j]) > 2.4)continue;
			if(obj1->VtxProb[j] < 0.001)continue;
			if (obj1->centrality/2 >= 0 && obj1->centrality/2 < cendiv){
				bool Acut = 1;
				
				if ((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j])) > 3.1415926){
					double A = 1-((2*3.1415926 - abs(obj1->PhiD1[j]) - abs(obj1->PhiD2[j]))/ 3.1415926);
					mc_reco_A_plot_24->Fill(A);
					if (A <= 0.001) Acut = 0;
				}
				else{
					double A = 1-((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j]))/ 3.1415926);
					mc_reco_A_plot_24->Fill(A);
					if (A <= 0.001) Acut = 0;
				}
				if (Acut) {
					mc_reco_Zmass_0_24->Fill(obj1->mass[j]);
				}
			}
			if (obj1->centrality/2 >= cendiv && obj1->centrality/2 < 100){
				bool Acut = 1;

				if ((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j])) > 3.1415926){
					double A = 1-((2*3.1415926 - abs(obj1->PhiD1[j]) - abs(obj1->PhiD2[j]))/ 3.1415926);
					mc_reco_A_plot_100->Fill(A);
					if (A <= 0.001) Acut = 0;
				}
				else{
					double A = 1-((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j]))/ 3.1415926);
					mc_reco_A_plot_100->Fill(A);
					if (A<=0.001) Acut = 0;
				}
				if (Acut){
					mc_reco_Zmass_24_100->Fill(obj1->mass[j]);
				}
			}
			//cout <<"PIDD1 is " << obj1->PIDD1[j] << " PIDD2 is " << obj1->PIDD2[j] << endl;
			// Why PIDD1 or PIDD2 == -77 ???
			//if (abs(obj1->PIDD1[j]) != 13 ||obj1->PIDD2[j] != 13) cout << "PID for reco D1 is " << obj1->PIDD1[j] << " PID for reco D2 is " << obj1->PIDD2[j] << endl;
			bool Acut = 1;
			if ((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j])) > 3.1415926){
				double A = 1-((2*3.1415926 - abs(obj1->PhiD1[j]) - abs(obj1->PhiD2[j]))/ 3.1415926);
				mc_reco_A_plot_all->Fill(A);
				if (A <= 0.001) Acut = 0;
			}
			else{
				double A = 1-((abs(obj1->PhiD1[j]) + abs(obj1->PhiD2[j]))/ 3.1415926);
				mc_reco_A_plot_all->Fill(A);
				if (A <= 0.001) Acut = 0;
			}
			if (Acut) mc_reco_Zmass->Fill(obj1->mass[j]);
		}
	}

	//Draw plot
	TFile *file = new TFile("./mc_bk_plot/W_histogram.root","UPDATE");
	file->cd();
	mc_reco_Zmass->Write("data_hist",2);
	mc_reco_Zmass_0_24->Write("data_hist_24",2);
	mc_reco_Zmass_24_100->Write("data_hist_100",2);
	mc_reco_A_plot_24->Write("data_A_24",2);
	mc_reco_A_plot_100->Write("data_A_100",2);
	mc_reco_A_plot_all->Write("data_A_ALL",2);
	file->Close();
	obj1->f1->Close();

}