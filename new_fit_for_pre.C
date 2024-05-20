#include "fitting_helper.h"
#include "MC_18.h"
void new_fit_for_pre(){
	fitting_helper *reco = new fitting_helper();
	fitting_helper *gen = new fitting_helper();
	cout <<"1" << endl;
	MC_18 *obj1 = new MC_18();

	reco->f1 = new TFile("./justplot.root","READ");
	gen->f1 = new TFile("./justplot.root","READ");

	cout << "2" << endl;

	reco->h_normalized_data[0] = (TH1D*) reco->f1->Get("mc_reco_hist");
	//gen->h_normalized_data[0] = (TH1D*) reco->f1->Get("mc_gen_hist");

		//Canvas
	TCanvas *c1 = new TCanvas("","",800,600);
	TCanvas *c2 = new TCanvas("","",800,600);

	//TPad
	TPad *p_1 = new TPad("p_1","Large",0,0.3,1,1);
	TPad *p_2 = new TPad("p_2","Small",0,0,1,0.3);
	TPad *p_3 = new TPad("p_3","Large",0,0.3,1,1);
	TPad *p_4 = new TPad("p_4","Small",0,0,1,0.3);

	TH1D *h_p_1 = new TH1D("h_p_1","Pull",120,60,120);

	obj1->Plot_and_Pull(reco->h_normalized_data[0],h_p_1,c1,p_1,p_2);

	c1->SaveAs("Temp.pdf");



	cout << "3" << endl;

	//reco->setupfittingfunction(0);
	//gen->setupfittingfunction(0);

	cout << "4" << endl;

	//reco->fitting(1);
	//gen->fitting(1);

	//reco->plotting(1,0,0);
	//gen->plotting(1,0,0);
}