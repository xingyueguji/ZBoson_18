#include "fitting_helper.h"
void get_fit(){
	//need 11 objs
	
	fitting_helper *obj[11];
	fitting_helper *obj1[11];
	fitting_helper *obj2[11];

	double Xcent[7] = {2.5,7.5,15,25,65,7.5,57.5}; 
	double xerr[7] = {0,0,0,0,0,0,0};
	double dMean[7] = {0,0,0,0,0,0,0};
	double yerr[7] = {0,0,0,0,0,0,0};
	for(int i=0; i<11; i++){
	obj[i] = new fitting_helper();
	obj1[i] = new fitting_helper();
	obj2[i] = new fitting_helper();
	}

	for (int i=0; i<7; i++){
		cout << "We are running centrality number " << i << endl;
		obj[i]->gethistograms();
		obj1[i]->gethistograms();
		obj2[i]->gethistograms();
		obj[i]->setupfittingfunction(i);
		obj1[i]->setupfittingfunction(i);
		obj2[i]->setupfittingfunction(i);
		//obj[i]->fitting(1);
		//obj1[i]->fitting(2);
		obj2[i]->fitting(3);
		//obj[i]->plotting(1,i);
		//obj1[i]->plotting(2,i);
		obj2[i]->plotting(3,i);
	}
	// obj[i] out of scope
	for (int i = 0; i < 7; i++){
		dMean[i] = obj2[i]->meanshift;
		yerr[i] = obj2[i]->meanshifterror;
		cout << "TEST!!!!!!" << obj2[i]->meanshift << endl;
	}
	TGraphErrors* g1 = new TGraphErrors(5,Xcent,dMean,xerr,yerr);

	TCanvas* c_2 = new TCanvas("","",800,600);
	c_2->cd();
	g1->Draw("AP");

	c_2->SaveAs("testdM.pdf");



}