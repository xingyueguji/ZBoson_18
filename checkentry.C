#include "fitting_helper.h"
void checkentry(bool rapiditycut = 1){
	fitting_helper *obj[11];
	fitting_helper *obj1[11];
	fitting_helper *obj2[11];
	for(int i=0; i<11; i++){
	obj[i] = new fitting_helper();
	obj[i]->gethistograms(rapiditycut);
	}

	TCanvas *c1 = new TCanvas("","",800,600);

	for(int i=0; i<11; i++){
		c1->cd();
		obj[i]->h_normalized_data[i]->Draw();
		if (!rapiditycut)c1->SaveAs(Form("./sancheck/data_hist_%i.pdf",i));
		if (rapiditycut)c1->SaveAs(Form("./sancheck/data_hist_%i_withycut.pdf",i));
		
	}

}