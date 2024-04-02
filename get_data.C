#include "MC_18.h"
void get_data(){
	TH1::SetDefaultSumw2();

	MC_18 *data = new MC_18();
	MC_18 *data_same_sign = new MC_18();

	data->SetupRootfile(4,0);
	data->SetupBranches(1);
	data_same_sign->SetupRootfile(5,0);
	data_same_sign->SetupBranches(1);


	TH1D *mass_array_data[11];
	TH1D *data_A[11];
	TH1D *mass_array_data_without_A_cut[11];

	TH1D *mass_array_data_same_sign[11];


	for (int i = 0; i < data->centarraysize; i++){
		mass_array_data[i] = new TH1D(Form("mass_array_data_%i",i),Form("mass_data_%i",data->cenhighlimit[i]),120,60,120);
		data_A[i] = new TH1D(Form("data_A_%i",i),Form("data_A_%i",data->cenhighlimit[i]),1000,0,1);
		mass_array_data_same_sign[i] = new TH1D(Form("mass_array_data_same_sign_%i",i),Form("mass_data_same_sign_%i",data_same_sign->cenhighlimit[i]),120,60,120);
		mass_array_data_without_A_cut[i] = new TH1D(Form("mass_array_data_noA_%i",i),Form("mass_data_noA_%i",data->cenhighlimit[i]),120,60,120);
	}

	int nentries = data->t1->GetEntries();
	cout << "We have " << nentries << " events in total" << endl;


	for (int i= 0; i<nentries; i++){
		if (i%100000 == 0) cout << "Event " << i << endl;
		//cout << "Event " << i << endl;
		data->t1->GetEntry(i);

		//Cut setup		

		for (int j=0;j<data->candSize; j++){
			if(data->mass[j] < data->masslowlimit || data->mass[j] > data->masshighlimit) continue;
			if(data->pTD1[j] < data->ptlowlimit) continue;
			if(data->pTD2[j] < data->ptlowlimit) continue;
			if(abs(data->EtaD1[j]) > 2.4)continue;
			if(abs(data->EtaD2[j]) > 2.4)continue;
			if(data->VtxProb[j] < 0.001)continue;
			int centbinpositioncounter[11] = {};
			data->CentBinSearching(centbinpositioncounter);
			for (int k=0; k<data->centarraysize; k++){
				if (centbinpositioncounter[k] != 0 ){
					bool Acut = 1;
					mass_array_data_without_A_cut[k]->Fill(data->mass[j]);
					if ((abs(data->PhiD1[j]) + abs(data->PhiD2[j])) > 3.1415926){
						double A = 1-((2*3.1415926 - abs(data->PhiD1[j]) - abs(data->PhiD2[j]))/ 3.1415926);
						data_A[k]->Fill(A);
						if (A <= 0.001) Acut = 0;
					}
					else{
						double A = 1-((abs(data->PhiD1[j]) + abs(data->PhiD2[j]))/ 3.1415926);
						data_A[k]->Fill(A);
						if (A <= 0.001) Acut = 0;
					}
					if (Acut) {
						mass_array_data[k]->Fill(data->mass[j]);
					}
				}
			}
		}
	}


	for (int i = 0; i < data_same_sign->t1->GetEntries(); i++){
		data_same_sign->t1->GetEntry(i);
		for (int j=0;j<data_same_sign->candSize; j++){
			if(data_same_sign->mass[j] < data_same_sign->masslowlimit || data_same_sign->mass[j] > data_same_sign->masshighlimit) continue;
			if(data_same_sign->pTD1[j] < data_same_sign->ptlowlimit) continue;
			if(data_same_sign->pTD2[j] < data_same_sign->ptlowlimit) continue;
			if(abs(data_same_sign->EtaD1[j]) > 2.4)continue;
			if(abs(data_same_sign->EtaD2[j]) > 2.4)continue;
			if(data_same_sign->VtxProb[j] < 0.001)continue;
			int centbinpositioncounter[11] = {};
			data_same_sign->CentBinSearching(centbinpositioncounter);
			for (int k=0; k<data_same_sign->centarraysize; k++){
				if (centbinpositioncounter[k] != 0 ){
						mass_array_data_same_sign[k]->Fill(data_same_sign->mass[j]);
					}
				}
			}
		}
	TFile *histogram_file = new TFile("./rootfile/data.root","UPDATE");
	histogram_file->cd();
	for (int i = 0; i < data->centarraysize; i++){
		mass_array_data[i]->Write("",2);
		data_A[i]->Write("",2);
		mass_array_data_same_sign[i]->Write("",2);
		mass_array_data_without_A_cut[i]->Write("",2);
	}
	histogram_file->Close();
	data->f1->Close();
}