#include "MC_18.h"

void MC_plot(bool isCut = 1){
	TH1::SetDefaultSumw2();
	MC_18 *obj1 = new MC_18();
	obj1->SetupRootfile();
	obj1->SetupBranches();

	//variable declaration
	bool eventnotpassed = 0; // 1 means not passed, 0 means passed cut

	//Histogram
	TH1D *mc_reco_Zmass = new TH1D("mc_reco_Zmass","mc_reco_Zmass",60,60,120);
	TH1D *h_p_1 = new TH1D("h_p_1","Pull",60,60,120);

	//Canvas
	TCanvas *c1 = new TCanvas("","",800,600);

	//TPad
	TPad *p_1 = new TPad("p_1","Large",0,0.3,1,1);
	TPad *p_2 = new TPad("p_2","Small",0,0,1,0.3);

	//Event loop

	int nentries = obj1->t1->GetEntries();
	cout << "We have " << nentries << " events in total" << endl;

	//ofstream out1("output.txt",std::ios::trunc);

	for (int i=0; i<nentries; i++){
		 //out1.seekp(0);
		 //out1 << "Updated content for iteration " << i << std::endl;
		if (i%1 == 0) cout << "Event " << i << endl;
		obj1->t1->GetEntry(i);
		cout << "TEST3" << endl;
		eventnotpassed = obj1->SetupCut();
		if (isCut && eventnotpassed) continue; // using cut and not passed, skip event

		for (int j=0;j<obj1->candSize; j++){
			cout << "TEST" << endl;
			if (obj1->mass[j] < 60 || obj1->mass[j] >= 120 ) continue;
			mc_reco_Zmass->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
		}

		//out1.flush();
		cout << "TEST2" << endl;
		
	}


	//Draw plot
	obj1->Plot_and_Pull(mc_reco_Zmass,h_p_1,c1,p_1,p_2);
	c1->SaveAs("./mcplot/mc_reco_Zmass.pdf");


}