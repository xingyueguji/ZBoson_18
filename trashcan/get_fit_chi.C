#include "fitting_helper.h"
void get_fit_chi(){
	//need 11 objs
	
	fitting_helper *obj[11];
	fitting_helper *obj1[11];
	fitting_helper *obj2[11];

	double Xcent[8] = {2.5, 7.5, 15,25,35,45,60,80}; 
	double xerr[8] = {0,0,0,0,0,0,0,0};
	double dMean[8] = {};
	double yerr[8] = {};

	for (int i=0; i<5; i++){
		cout << "We are running centrality number " << i << endl;
		obj[i] = new fitting_helper();
		obj1[i] = new fitting_helper();
		obj2[i] = new fitting_helper();
		obj[i]->gethistograms();
		obj1[i]->gethistograms();
		obj2[i]->gethistograms();
		obj[i]->setupfittingfunction(i);
		obj1[i]->setupfittingfunction(i);
		obj2[i]->setupfittingfunction(i);
		obj[i]->fitting(1);
		//obj1[i]->fitting(2);
		//obj2[i]->fitting(3);
		obj[i]->plotting(1,i);
		//obj1[i]->plotting(2,i);
		//obj2[i]->plotting(3,i);
	}
	for (int i = 0; i < 8; i++){
		dMean[i] = obj[i] -> meanarray[i];
		yerr[i] = obj[i] -> meanerrorarray[i];
	}
	TGraphErrors* g1 = new TGraphErrors(8,Xcent,dMean,xerr,yerr);

	TCanvas* c_2 = new TCanvas("","",800,600);
	g1->Draw("AP");

	c_2->SaveAs("testdM.pdf");



}