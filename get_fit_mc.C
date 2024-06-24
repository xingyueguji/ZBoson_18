#include "fitting_helper.h"
void get_fit_mc(int type = 3, bool iseff = 1, bool w2unbin = 1){
	//need 11 objs
	// type 1 = raw + reco
	// type 2 = ycut + reco
	// type 3 = etacut + reco
	// type 4 = raw + gen
	// type 5 = ycut + gen
	// type 6 = etacut + gen

	// iseff only affect reco, not gen!
	
	fitting_helper *obj[11];

	double Xcent[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr[7] = {0,0,0,0,0,0,0};
	double Ycent[7] = {0,0,0,0,0,0,0};
	double yerr[7] = {0,0,0,0,0,0,0};

	double Xcent1[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr1[7] = {0,0,0,0,0,0,0};
	double Ycent1[7] = {0,0,0,0,0,0,0};
	double yerr1[7] = {0,0,0,0,0,0,0};

	double Xcent2[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr2[7] = {0,0,0,0,0,0,0};
	double Ycent2[7] = {0,0,0,0,0,0,0};
	double yerr2[7] = {0,0,0,0,0,0,0};

	double Xcent3[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr3[7] = {0,0,0,0,0,0,0};
	double Ycent3[7] = {0,0,0,0,0,0,0};
	double yerr3[7] = {0,0,0,0,0,0,0};

	double Xcent4[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr4[7] = {0,0,0,0,0,0,0};
	double Ycent4[7] = {0,0,0,0,0,0,0};
	double yerr4[7] = {0,0,0,0,0,0,0};

	for(int i=0; i<11; i++){
	obj[i] = new fitting_helper();
	/*if (type == 3 && i == 0){
		obj[i]->bwmeanrecoip = 91.1373;
		obj[i]->widthrecoip = 2.5159;
		obj[i]->sigmarecoip = 0.7414;
		obj[i]->alpharecoip = 1.8067;
		obj[i]->nrecoip = 0.9393;

	}
	if (type == 3 && i == 1){
		obj[i]->bwmeanrecoip = 91.0799;
		obj[i]->widthrecoip = 2.4176;
		obj[i]->sigmarecoip = 0.8592;
		obj[i]->alpharecoip = 1.7849;
		obj[i]->nrecoip = 0.9414;
		
	}
	if (type == 3 && i == 2){
		obj[i]->bwmeanrecoip = 91.1429;
		obj[i]->widthrecoip = 2.4646;
		obj[i]->sigmarecoip = 0.7475;
		obj[i]->alpharecoip = 1.7618;
		obj[i]->nrecoip = 0.9773;
		
	}
	if (type == 3 && i == 3){
		obj[i]->bwmeanrecoip = 91.0999;
		obj[i]->widthrecoip = 2.5070;
		obj[i]->sigmarecoip = 0.8213;
		obj[i]->alpharecoip = 1.8427;
		obj[i]->nrecoip = 0.9274;
		
	}
	if (type == 3 && i == 4){
		obj[i]->bwmeanrecoip = 91.0976;
		obj[i]->widthrecoip = 2.4336;
		obj[i]->sigmarecoip = 0.8589;
		obj[i]->alpharecoip = 1.8165;
		obj[i]->nrecoip = 0.9257;
		
	}
	if (type == 3 && i == 5){
		obj[i]->bwmeanrecoip = 91.1218;
		obj[i]->widthrecoip = 2.4730;
		obj[i]->sigmarecoip = 0.7834;
		obj[i]->alpharecoip = 1.7989;
		obj[i]->nrecoip = 0.9377;
		
	}
	if (type == 3 && i == 6){
		obj[i]->bwmeanrecoip = 91.1087;
		obj[i]->widthrecoip = 2.4791;
		obj[i]->sigmarecoip = 0.8039;
		obj[i]->alpharecoip = 1.8041;
		obj[i]->nrecoip = 0.9473;
		
	}*/
	}

	for (int i=0; i<7; i++){
		cout << "We are running centrality number " << i << endl;
		obj[i]->getmcdataset();

		obj[i]->setupfittingfunctionmc();

		obj[i]->fittingmc(type,iseff,w2unbin,i);

		obj[i]->plottingmc(type,i,iseff,w2unbin);

	}

	for (int i=0; i<7; i++){
		Ycent[i] = obj[i]->meanshift;
		yerr[i] = obj[i]->meanshifterror;

		Ycent1[i] = obj[i]->widthshift;
		yerr1[i] = obj[i]->widtherror;

		if (type < 4){
			Ycent2[i] = obj[i]-> cbsigma->getVal();
			yerr2[i] = obj[i]-> cbsigma->getError();

			Ycent3[i] = obj[i] -> cbalpha->getVal();
			yerr3[i] = obj[i] -> cbalpha->getError();

			Ycent4[i] = obj[i]-> cbn->getVal();
			yerr4[i] = obj[i]-> cbn->getError();
		}
		if (type >= 4){
			Ycent2[i] = obj[i]-> cbsigmagen->getVal();
			yerr2[i] = obj[i]-> cbsigmagen->getError();

			Ycent3[i] = obj[i] -> cbalphagen->getVal();
			yerr3[i] = obj[i] -> cbalphagen->getError();

			Ycent4[i] = obj[i]-> cbngen->getVal();
			yerr4[i] = obj[i]-> cbngen->getError();
		}


	}

	TGraphErrors* g1 = new TGraphErrors(7,Xcent,Ycent,xerr,yerr);
	TGraphErrors* g2 = new TGraphErrors(7,Xcent1,Ycent1,xerr1,yerr1);
	TGraphErrors* g3 = new TGraphErrors(7,Xcent2,Ycent2,xerr2,yerr2);
	TGraphErrors* g4 = new TGraphErrors(7,Xcent3,Ycent3,xerr,yerr3);
	TGraphErrors* g5 = new TGraphErrors(7,Xcent4,Ycent4,xerr1,yerr4);

	obj[1] -> plotTGraphmc(g1,type,iseff,w2unbin,1,1);
	obj[1] -> plotTGraphmc(g2,type,iseff,w2unbin,2,1);
	obj[1] -> plotTGraphmc(g3,type,iseff,w2unbin,3,1);
	obj[1] -> plotTGraphmc(g4,type,iseff,w2unbin,4,1);
	obj[1] -> plotTGraphmc(g5,type,iseff,w2unbin,5,1);

}