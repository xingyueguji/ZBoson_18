#include "chisquaretest.h"

void newchisquaretest(int type = 1, bool isbk = 0){

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	TString mcfile;
	TString datafile;
	TString bkfile;

	if (type == 0) bkfile = "./rootfile/normalized/efffile.root";
	if (type == 1) bkfile = "./rootfile/normalized/etacut_eff_file.root";

	datafile = "./rootfile/data_file.root";
	mcfile = "./rootfile/modified_signal.root";

	chisquaretest* ovo = new chisquaretest(mcfile,datafile,bkfile,type);



	//ovo->bincontentcheck(isbk);
	ovo->RebinAll(4);
	ovo->calculatechisq(isbk);
	ovo->plottingandformatting(type,isbk);

}