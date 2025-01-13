#include "../../pp_18/chisquaretest.h"

void testtemplate()
{
    chisquaretest* ovo = new chisquaretest();

    TFile *f1 = new TFile("new_template_pp_reco_gen_raw_test_-1_1_-1_1.root", "READ");
    TFile *f2 = new TFile("new_template_pp_reco_gen_raw_test_-1_1_0_2.root","READ");

    TFile *pp = new TFile("../../pp_18/new_pp_data_file_stability_readonly.root","READ");

    TH1D* template1;
    TH1D* template2;
    TH1D* h_pp;

    template1 = (TH1D*) f1->Get("mass_array_with_eff_template_10_10_10");
    template2 = (TH1D*) f2->Get("mass_array_with_eff_template_10_0_10");
    h_pp = (TH1D*) pp->Get("h_mass_array_raw_0");

    ovo->areanormalize(template1);
    ovo->areanormalize(template2);
    ovo->areanormalize(h_pp);
    cout << "first chi2 is" << ovo->myownfunctionchi2(h_pp, template1) << " second chi2 is " << ovo->myownfunctionchi2(h_pp,template2) << endl;



    /*template1->Add(template2,-1);

    TCanvas* c1 = new TCanvas("","",800,600);

    template1->Draw();

    c1->SaveAs("test.png");*/

}