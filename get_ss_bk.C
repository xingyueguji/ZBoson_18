void get_ss_bk(){
	TH1::SetDefaultSumw2();
	int nbins = 120;
	TFile *mc_file = new TFile("./rootfile/mc_bk.root","UPDATE");
	TFile *data_file = new TFile("./rootfile/data.root","UPDATE");

	TH1D* data_os[11];
	TH1D* mc_os[11];
	TH1D* data_ss[11];

	TH1D* same_sign_ratio[11];

	for (int i = 0; i < 11; i++){
		data_os[i] = (TH1D*) data_file->Get(Form("mass_array_data_noA_%i",i));
		mc_os[i] = (TH1D*) mc_file->Get(Form("mass_array_signal_%i",i));
		data_ss[i] = (TH1D*) data_file->Get(Form("mass_array_data_same_sign_%i",i));
		same_sign_ratio[i] = new TH1D(Form("same_sign_ratio_%i",i),Form("same_sign_ratio_%i",i),120,60,120);
	}
	for (int i = 0; i < 11; i++){
		for (int j = 1; j<=120; j++){
			double content_mc_os = mc_os[i]->GetBinContent(j);
			double content_data_ss = data_ss[i]->GetBinContent(j);
			double content_data_os = data_os[i]->GetBinContent(j);
			double ratio;
			if (content_data_os >0 ) ratio = content_data_ss / content_data_os;
			else ratio = 0;
			//if (i == 7) cout << "Ratio is " <<  ratio << endl;
			double estimate_mc_ss = content_mc_os * ratio;
			if (estimate_mc_ss >= 0 )same_sign_ratio[i]->SetBinContent(j,estimate_mc_ss);
			//if (i == 7) cout << "Bin content is " <<  estimate_mc_ss << endl;
		}
	}
	mc_file->cd();
	for (int i=0; i<11; i++)same_sign_ratio[i]->Write(Form("mc_estimate_ss_%i",i),2);
	mc_file->Close();
	data_file->Close();

}