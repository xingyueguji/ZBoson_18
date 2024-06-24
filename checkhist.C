void checkhist(){
	TString raw_file = "./rootfile/normalized/rawfile.root";
	TString ycut_file = "./rootfile/normalized/ycutfile.root";
	TString eff_file = "./rootfile/normalized/efffile.root";
	TString ycut_eff_file = "./rootfile/normalized/ycut_eff_file.root";
	TString etacut_file = "./rootfile/normalized/etacut_file.root";
	TString etacut_eff_file = "./rootfile/normalized/etacut_eff_file.root";

	TString data_file = "./rootfile/data_file.root";
	TString mc_file = "./rootfile/mc_signal.root";

	TFile *f1;
	TFile *f2;

	f1 = new TFile(etacut_eff_file,"READ");
	f2 = new TFile(mc_file,"READ");

	TH1D *hist;
	RooDataSet *set;

	for (int i=0; i<7; i++){
		hist = (TH1D*)f2->Get(Form("mass_array_witheta_witheff_%i",i));
		set = (RooDataSet*)f2->Get(Form("rooreco_eta_eff_%i",i));

		RooRealVar *x = new RooRealVar("roomass","roomass",60,120);
		RooPlot* frame = x->frame(RooFit::Title("Z mass fit"));

		set->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(120),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));

		TCanvas *c1 = new TCanvas("","",800,400);

		c1->Divide(2,1);
		c1->cd(1);

		hist->Draw();

		c1->cd(2);

		frame->Draw();

		c1->SaveAs(Form("./verifyhistogram/testhist_%i.pdf",i));
	}



}