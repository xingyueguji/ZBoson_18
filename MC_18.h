#include "TMath.h"
class MC_18{
	public:

	MC_18();
	~MC_18();
	
	void SetupRootfile();
	void SetupBranches();
	bool SetupCut();
	void Plot_and_Pull(TH1D* hist, TH1D* hist2, TCanvas* c1, TPad* p1, TPad* p2);
	void Plot_and_Pull_gen(TH1D* hist, TH1D* hist2, TCanvas* c1, TPad* p1, TPad* p2);
    void Calc_Z_gen(TH1D* h1);
	

	//TString MC_File_Path = "/eos/cms/store/group/phys_heavyions/abaty/Zmumu_Samples_MC/VertexCompositeTree_DYJetsToLL_TuneCP5_HydjetDrumMB_Pythia8_HINPbPbAutumn18_DiMuMC_20190808.root";
	TString MC_File_Path = "./VertexCompositeTree_DYJetsToLL_TuneCP5_HydjetDrumMB_Pythia8_HINPbPbAutumn18_DiMuMC_20190808.root";
	TFile *f1;
	TDirectoryFile *dir1;
	TTree *t1;


    //These are reco level;
	UInt_t candSize = 0;
	Short_t centrality = 0;
	Float_t pT[5000] = {};
	Float_t eta[5000] = {};
	Float_t phi[5000] = {};
	Float_t mass[5000] = {};
	Float_t VtxProb[5000] = {};
	Float_t pTD1[5000] = {};
	Float_t EtaD1[5000] = {};
	Float_t PhiD1[5000] = {};
	Short_t chargeD1[5000] = {};
	Float_t pTD2[5000] = {};
	Float_t EtaD2[5000] = {};
	Float_t PhiD2[5000] = {};
	Short_t chargeD2[5000] = {};
	vector<float> *weightLHE_gen = 0;

	//These are I believe gen level
	UInt_t candSize_gen = 0;
	Float_t weight_gen = 0;
	Short_t DecayID_gen[5000] = {};
	Float_t pTD1_gen[5000] = {};
	Float_t EtaD1_gen[5000] = {};
	Float_t PhiD1_gen[5000] = {};
	Float_t pTD2_gen[5000] = {};
	Float_t EtaD2_gen[5000] = {};
	Float_t PhiD2_gen[5000] = {};
	Int_t PIDD2_gen[5000] = {};
	Int_t PIDD1_gen[5000] = {};


	Int_t nentries = 0;

	


};

MC_18::MC_18(){

}

MC_18::~MC_18(){
	f1->Close();
	delete f1;
	delete t1;
	delete dir1;
}

void MC_18::SetupRootfile(){
	f1 = new TFile(MC_File_Path,"READ");
	if (f1->IsOpen()!=1){
		cout << "ROOTfile " << MC_File_Path << " Cannot be opened! " << endl; 
	}
	dir1 = (TDirectoryFile*)f1->Get("dimucontana_mc;1");
	if (!dir1){
		cout << "TDirectory dimucontana_mc cannot be opened! " << endl; 
	}
	t1 = (TTree*) dir1->Get("VertexCompositeNtuple;1");
	if (!t1){
		cout << "TTree VertexCompositeNtuple;1 cannot be opened! " << endl;
	}
}
void MC_18::SetupBranches(){
	//Enable required branches
	t1->SetBranchStatus("*",0);
	t1->SetBranchStatus("candSize",1);
	t1->SetBranchStatus("centrality",1);
	t1->SetBranchStatus("pT",1);
	t1->SetBranchStatus("eta",1);
	t1->SetBranchStatus("phi",1);
	t1->SetBranchStatus("mass",1);
	t1->SetBranchStatus("VtxProb",1);
	t1->SetBranchStatus("pTD1",1);
	t1->SetBranchStatus("EtaD1",1);
	t1->SetBranchStatus("PhiD1",1);
	t1->SetBranchStatus("chargeD1",1);
	t1->SetBranchStatus("pTD2",1);
	t1->SetBranchStatus("EtaD2",1);
	t1->SetBranchStatus("PhiD2",1);
	t1->SetBranchStatus("chargeD2",1);
	t1->SetBranchStatus("weightLHE_gen",1);
	t1->SetBranchStatus("candSize_gen",1);
	t1->SetBranchStatus("weight_gen",1);
	t1->SetBranchStatus("DecayID_gen",1);
	t1->SetBranchStatus("pTD1_gen",1);
	t1->SetBranchStatus("EtaD1_gen",1);
	t1->SetBranchStatus("PhiD1_gen",1);
	t1->SetBranchStatus("pTD2_gen",1);
	t1->SetBranchStatus("EtaD2_gen",1);
	t1->SetBranchStatus("PhiD2_gen",1);
	t1->SetBranchStatus("PIDD2_gen",1);
	t1->SetBranchStatus("PIDD1_gen",1);
	//Link to variables
	t1->SetBranchAddress("candSize",&candSize);
	t1->SetBranchAddress("centrality",&centrality);
	t1->SetBranchAddress("pT",pT);
	t1->SetBranchAddress("eta",eta);
	t1->SetBranchAddress("phi",phi);
	t1->SetBranchAddress("mass",mass);
	t1->SetBranchAddress("VtxProb",VtxProb);
	t1->SetBranchAddress("pTD1",pTD1);
	t1->SetBranchAddress("EtaD1",EtaD1);
	t1->SetBranchAddress("PhiD1",PhiD1);
	t1->SetBranchAddress("chargeD1",chargeD1);
	t1->SetBranchAddress("pTD2",pTD2);
	t1->SetBranchAddress("EtaD2",EtaD2);
	t1->SetBranchAddress("PhiD2",PhiD2);
	t1->SetBranchAddress("chargeD2",chargeD2);
	t1->SetBranchAddress("weightLHE_gen",&weightLHE_gen);
	t1->SetBranchAddress("candSize_gen",&candSize_gen);
	t1->SetBranchAddress("weight_gen",&weight_gen);
	t1->SetBranchAddress("DecayID_gen",DecayID_gen);
	t1->SetBranchAddress("pTD1_gen",pTD1_gen);
	t1->SetBranchAddress("EtaD1_gen",EtaD1_gen);
	t1->SetBranchAddress("PhiD1_gen",PhiD1_gen);
	t1->SetBranchAddress("pTD2_gen",pTD2_gen);
	t1->SetBranchAddress("EtaD2_gen",EtaD2_gen);
	t1->SetBranchAddress("PhiD2_gen",PhiD2_gen);
	t1->SetBranchAddress("PIDD2_gen",PIDD2_gen);
	t1->SetBranchAddress("PIDD1_gen",PIDD1_gen);
}
bool MC_18::SetupCut(){
	int isPassed = 0;
	for (int i = 0; i < candSize; i++){
		if (pTD1[i] > 20 && TMath::Abs(EtaD1[i]) < 2.4 && pTD2[i] > 20 && TMath::Abs(EtaD2[i]) < 2.4 && VtxProb[i] > 0.001) isPassed++;
	}

	if (isPassed == candSize){
		return false;
	}
	else{
		return true;
	}
}
void MC_18::Plot_and_Pull(TH1D* hist, TH1D* hist2, TCanvas* c1, TPad* p1, TPad* p2){
	gStyle->SetOptFit(1);
	hist2->GetXaxis()->SetTitle("Invariant Mass(GeV)");
	c1->cd();
	TF1 *fit = new TF1("voigt","[0]*TMath::Voigt(x-[1]-91.1876,[2],[3]+2.4955)",60,120);
	int numBins = hist->GetNbinsX();
	int maxBin = hist->GetMaximumBin();
	fit->SetParameter(0,maxBin);
	//fit->SetParameter(1,h_m_1->GetMean()-91.1876);
	fit->SetParameter(1,0);
	fit->SetParameter(2,hist->GetStdDev());
	fit->SetLineWidth(1);
	fit->SetParameter(3,2.4955);
	p1->SetMargin(0.1,0.1,0,0.1);
	p2->SetMargin(0.1,0.1,0.2,0);
	p1->Draw();
	p2->Draw();
	p1->cd();
	//p1->SetLogy();
	hist->Fit(fit,"LQR");
	hist->Draw();
	Double_t integral = hist->Integral(1,numBins);
	TLatex *text = new TLatex();
	text->SetNDC();
	text->DrawLatex(0.2,0.8,Form("# of entries is %i",int(integral)));
	p2->cd();
	for (int i=1; i<=numBins; i++){
		double bincenter = hist->GetXaxis()->GetBinCenter(i);
		double fitresult = fit->Eval(bincenter);
		double binerror = hist->GetBinError(i);
		double bincontent = hist->GetBinContent(i);
		double pullvalue = (bincontent - fitresult) / binerror;
		hist2->SetBinContent(i,pullvalue);
		hist2->SetBinError(i,1);
	}
	hist2->GetYaxis()->SetLabelSize(0.08);
	hist2->GetXaxis()->SetLabelSize(0.1);
	hist2->GetXaxis()->SetTitleSize(0.1);
	hist2->SetStats(0);
	hist2->SetTitle("");
	hist2->GetYaxis()->SetTitle("Pull");
	hist2->GetYaxis()->SetTitleSize(0.15);
	hist2->GetYaxis()->SetTitleOffset(0.2);
	hist2->GetXaxis()->SetTitleOffset(0.9);
	hist2->GetXaxis()->SetTitle("Invariant Mass(GeV)");
	hist2->GetYaxis()->CenterTitle();
	hist2->Draw("P HIST E1");

}
void MC_18::Plot_and_Pull_gen(TH1D* hist, TH1D* hist2, TCanvas* c1, TPad* p1, TPad* p2){
	gStyle->SetOptFit(1);
	hist2->GetXaxis()->SetTitle("Invariant Mass(GeV)");
	c1->cd();
	TF1 *fit = new TF1("voigt","[0]*TMath::Voigt(x-[1]-91.1876,[2],[3]+2.4955)",61.1876,121.1876);
	//TF1 *fit = new TF1("BreitWigner","[0]*TMath::BreitWigner(x, [1]+91.1876, [2]+2.4955)",60,120);
	int numBins = hist->GetNbinsX();
	int maxBin = hist->GetMaximumBin();
	fit->SetParameter(0,20000);
	//fit->SetParameter(1,h_m_1->GetMean()-91.1876);
	fit->SetParameter(1,0);
	fit->SetLineWidth(1);
	fit->SetParameter(2,0);
	p1->SetMargin(0.1,0.1,0,0.1);
	p2->SetMargin(0.1,0.1,0.2,0);
	p1->Draw();
	p2->Draw();
	p1->cd();
	//p1->SetLogy();
	hist->Fit(fit,"LQR");
	hist->Draw();
	Double_t integral = hist->Integral(1,numBins);
	TLatex *text = new TLatex();
	text->SetNDC();
	text->DrawLatex(0.2,0.8,Form("# of entries is %i",int(integral)));
	p2->cd();
	for (int i=1; i<=numBins; i++){
		double bincenter = hist->GetXaxis()->GetBinCenter(i);
		double fitresult = fit->Eval(bincenter);
		double binerror = hist->GetBinError(i);
		double bincontent = hist->GetBinContent(i);
		double pullvalue = (bincontent - fitresult) / binerror;
		hist2->SetBinContent(i,pullvalue);
		hist2->SetBinError(i,1);
	}
	hist2->GetYaxis()->SetLabelSize(0.08);
	hist2->GetXaxis()->SetLabelSize(0.1);
	hist2->GetXaxis()->SetTitleSize(0.1);
	hist2->SetStats(0);
	hist2->SetTitle("");
	hist2->GetYaxis()->SetTitle("Pull");
	hist2->GetYaxis()->SetTitleSize(0.15);
	hist2->GetYaxis()->SetTitleOffset(0.2);
	hist2->GetXaxis()->SetTitleOffset(0.9);
	hist2->GetXaxis()->SetTitle("Invariant Mass(GeV)");
	hist2->GetYaxis()->CenterTitle();
	hist2->Draw("P HIST E1");

}
void MC_18::Calc_Z_gen(TH1D* h1){
	for (int i=0; i<this->candSize_gen; i++){
		TLorentzVector muon1;
		if (abs(PIDD2_gen[i]) != 13 || abs(PIDD1_gen[i]) !=13) continue;
		Float_t pt_muon_1 = this->pTD1_gen[i];
		Float_t pt_muon_2 = this->pTD2_gen[i];
		Float_t eta_muon_1 = this->EtaD1_gen[i];
		Float_t eta_muon_2 = this->EtaD2_gen[i];
		if (pt_muon_1 < 20 || pt_muon_2 < 20 || abs(eta_muon_1) > 2.4 || abs(eta_muon_2) > 2.4) continue;
		if (this->DecayID_gen[i] != 15 && this->DecayID_gen[i] != 23) cout << "??????????????" << this->DecayID_gen[i] << endl;
		if (this->DecayID_gen[i] != 23) continue;

		Float_t phi_muon_1 = this->PhiD1_gen[i];
		Float_t phi_muon_2 = this->PhiD2_gen[i];
		muon1.SetPtEtaPhiM(pt_muon_1,eta_muon_1,phi_muon_1,0.105);
		TLorentzVector muon2;
		muon2.SetPtEtaPhiM(pt_muon_2,eta_muon_2,phi_muon_2,0.105);
		TLorentzVector ZBoson = muon1 + muon2;
		Float_t M = ZBoson.M();
		if (M < 61.1876 || M >= 121.1876) continue;
		h1->Fill(M,this->weightLHE_gen->at(1080)/10000);
	}
}
