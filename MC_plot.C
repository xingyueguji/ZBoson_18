#include "MC_18.h"

void MC_plot(bool isCut = 1){
	TH1::SetDefaultSumw2();
	MC_18 *obj1 = new MC_18();
	obj1->SetupRootfile();
	obj1->SetupBranches();

	//variable declaration
	bool eventnotpassed = 0; // 1 means not passed, 0 means passed cut

	//Histogram
	TH1D *mc_reco_Zmass = new TH1D("mc_reco_Zmass","mc_reco_Zmass",120,60,120);
	TH1D *h_p_1 = new TH1D("h_p_1","Pull",120,60,120);
	TH1D *mc_gen_Zmass = new TH1D("mc_gen_Zmass","mc_gen_Zmass",120,60,120);
	TH1D *h_p_2 = new TH1D("h_p_2","Pull",120,60,120);

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

	for (int i=0; i<nentries; i++){
		if (i%100000 == 0) cout << "Event " << i << endl;
		obj1->t1->GetEntry(i);
		eventnotpassed = obj1->SetupCut();
		if (isCut && eventnotpassed) continue; // using cut and not passed, skip event

		for (int j=0;j<obj1->candSize; j++){
			if (obj1->mass[j] < 60 || obj1->mass[j] >= 120 || obj1->DecayID_gen[j] == 15) continue; 
			mc_reco_Zmass->Fill(obj1->mass[j],obj1->weightLHE_gen->at(1080)/10000.0);
		}
		obj1->Calc_Z_gen(mc_gen_Zmass);
	}


	//Draw plot
	obj1->Plot_and_Pull(mc_reco_Zmass,h_p_1,c1,p_1,p_2);
	obj1->Plot_and_Pull_gen(mc_gen_Zmass,h_p_2,c2,p_3,p_4);
	c1->SaveAs("./mcplot/mc_reco_Zmass.pdf");
	c2->SaveAs("./mcplot/mc_gen_Zmass.pdf");


}