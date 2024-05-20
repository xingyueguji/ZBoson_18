#include "fitting_helper.h"
void get_fit(int opt = 1){
	//need 11 objs
	// opt 1 == raw
	// opt 2 == y
	// opt 3 == eff
	// opt 4 == y+eff
	// opt 5 == eta
	// opt 6 == eta+eff
	
	fitting_helper *obj[11];
	fitting_helper *obj1[11];
	fitting_helper *obj2[11];

	double Xcent[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr[7] = {0,0,0,0,0,0,0};
	double dMean[7] = {0,0,0,0,0,0,0};
	double yerr[7] = {0,0,0,0,0,0,0};

	double Xcent1[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr1[7] = {0,0,0,0,0,0,0};
	double dMean1[7] = {0,0,0,0,0,0,0};
	double yerr1[7] = {0,0,0,0,0,0,0};

	double Xcent2[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerr2[7] = {0,0,0,0,0,0,0};
	double dMean2[7] = {0,0,0,0,0,0,0};
	double yerr2[7] = {0,0,0,0,0,0,0};

	double Xcentdsig[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerrdsig[7] = {0,0,0,0,0,0,0};
	double dMeandsig[7] = {0,0,0,0,0,0,0};
	double yerrdsig[7] = {0,0,0,0,0,0,0};

	double Xcentdsig1[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerrdsig1[7] = {0,0,0,0,0,0,0};
	double dMeandsig1[7] = {0,0,0,0,0,0,0};
	double yerrdsig1[7] = {0,0,0,0,0,0,0};

	double Xcentdsig2[7] = {5,15,25,40,75,7.5,57.5}; 
	double xerrdsig2[7] = {0,0,0,0,0,0,0};
	double dMeandsig2[7] = {0,0,0,0,0,0,0};
	double yerrdsig2[7] = {0,0,0,0,0,0,0};

	for(int i=0; i<11; i++){
	obj[i] = new fitting_helper();
	obj1[i] = new fitting_helper();
	obj2[i] = new fitting_helper();
	}

	for (int i=0; i<7; i++){
		cout << "We are running centrality number " << i << endl;
		obj[i]->gethistograms(opt);
		obj1[i]->gethistograms(opt);
		obj2[i]->gethistograms(opt);
		obj[i]->setupfittingfunction(i);
		obj1[i]->setupfittingfunction(i);
		obj2[i]->setupfittingfunction(i);
		obj[i]->fitting(1);
		obj1[i]->fitting(2);
		obj2[i]->fitting(3);
		obj[i]->plotting(opt,1,i);
		obj1[i]->plotting(opt,2,i);
		obj2[i]->plotting(opt,3,i);
	}
	// obj[i] out of scope
	for (int i = 0; i < 7; i++){
		dMean[i] = obj[i]->meanshift;
		yerr[i] = obj[i]->meanshifterror;
		dMean1[i] = obj1[i]->meanshift;
		yerr1[i] = obj1[i]->meanshifterror;
		dMean2[i] = obj2[i]->meanshift;
		yerr2[i] = obj2[i]->meanshifterror;

		dMeandsig[i] = obj[i]->widthshift;
		yerrdsig[i] = obj[i]->widtherror;
		dMeandsig1[i] = obj1[i]->widthshift;
		yerrdsig1[i] = obj1[i]->widtherror;
		dMeandsig2[i] = obj2[i]->widthshift;
		yerrdsig2[i] = obj2[i]->widtherror;
		
	}
	TGraphErrors* g1 = new TGraphErrors(7,Xcent,dMean,xerr,yerr);
	TGraphErrors* g2 = new TGraphErrors(7,Xcent1,dMean1,xerr1,yerr1);
	TGraphErrors* g3 = new TGraphErrors(7,Xcent2,dMean2,xerr2,yerr2);

	TGraphErrors* g1sig = new TGraphErrors(7,Xcentdsig,dMeandsig,xerrdsig,yerrdsig);
	TGraphErrors* g2sig = new TGraphErrors(7,Xcentdsig1,dMeandsig1,xerrdsig1,yerrdsig1);
	TGraphErrors* g3sig = new TGraphErrors(7,Xcentdsig2,dMeandsig2,xerrdsig2,yerrdsig2);

	obj[1]->plotTGraph(g1,g2,g3,1,opt);
	obj[1]->plotTGraph(g1sig,g2sig,g3sig,2,opt);

}