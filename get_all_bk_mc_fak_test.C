#include "MC_18.h"
#include "MuonTnP.h"
#include "MCReweight.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"
#include "TEfficiency.h"
#include "ptReweightSpectrum.h"
#include "plotting_helper.h"

// C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

void plotandsave(TCanvas *c1, TH1D *h1, std::string opt, const char *savepath, const char *addtionaltext)
{
    c1->cd();
    if (opt == "eta" || opt == "phi")
    {
        c1->SetLogy(true);
    }
    else
    {
        c1->SetLogy(false);
    }

    h1->GetYaxis()->SetTitle("Events");
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleSize(0.05);

    if (opt == "mass")
    {
        h1->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
    }
    if (opt == "pT")
    {
        h1->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)");
    }
    if (opt == "eta")
    {
        h1->GetXaxis()->SetTitle("#eta^{#mu}");
    }
    if (opt == "phi")
    {
        h1->GetXaxis()->SetTitle("#phi^{#mu}");
    }

    h1->SetMarkerSize(0.5);
    h1->Draw("");

    TPaveText *pave = new TPaveText(0.65, 0.75, 0.90, 0.88, "NDC");
    pave->SetBorderSize(0);
    pave->SetFillColor(0);
    pave->SetTextFont(42); // CMS font
    pave->SetTextSize(0.04);
    pave->SetTextAlign(12); // Align left

    // Add text lines
    pave->AddText(addtionaltext);
    pave->Draw("SAME");

    c1->SaveAs(savepath);
    c1->Clear();
}

bool checkcut(int j, MC_18 *s)
{
    if (s->mass[j] < s->masslowlimit || s->mass[j] > s->masshighlimit)
        return false;
    if (s->pTD1[j] < s->ptlowlimit)
        return false;
    if (s->pTD2[j] < s->ptlowlimit)
        return false;
    if (abs(s->EtaD1[j]) > 2.4)
        return false;
    if (abs(s->EtaD2[j]) > 2.4)
        return false;
    if (s->VtxProb[j] < 0.001)
        return false;
    if (abs(s->y[j]) > 2.4)
        return false;
    if (!(s->tightMuon1Cut(j, "POG") && s->tightMuon2Cut(j, "POG")))
        return false;

    bool isDaughter1Trigger = s->trigMuon1->at(6).at(j); // 6 is HLT_HIL3Mu12;
    bool isDaughter2Trigger = s->trigMuon2->at(6).at(j);
    if (!(isDaughter1Trigger || isDaughter2Trigger))
        return false;

    bool isOppositeSign = s->chargeD1[j] != s->chargeD2[j];
    if (!isOppositeSign)
    {
        return false;
    }
    return true;
}

bool checkduplicate(int j, MC_18 *s)
{
    for (int i = 0; i < s->candSize; i++)
    {
        if (i == j)
            continue;
        if (s->pTD1[j] == s->pTD1[i])
        {
            return true;
        }
        else if (s->pTD1[j] == s->pTD2[i])
        {
            return true;
        }
    }
    return false;
}

Double_t outputduplicate(int j, MC_18 *s, std::string type)
{
    for (int i = 0; i < s->candSize; i++)
    {
        if (i == j)
            continue;

        if (type == "pT")
        {
            if (s->pTD1[j] == s->pTD1[i])
            {
                return s->pTD1[j];
            }
            else if (s->pTD1[j] == s->pTD2[i])
            {
                return s->pTD1[j];
            }
        }

        if (type == "eta")
        {
            if (s->EtaD1[j] == s->EtaD1[i])
            {
                return s->EtaD1[j];
            }
            else if (s->EtaD1[j] == s->EtaD2[i])
            {
                return s->EtaD1[j];
            }
        }

        if (type == "phi")
        {
            if (s->PhiD1[j] == s->PhiD1[i])
            {
                return s->PhiD1[j];
            }
            else if (s->PhiD1[j] == s->PhiD2[i])
            {
                return s->PhiD1[j];
            }
        }
    }
    return -99;
}

void get_all_bk_mc_fak_test(int type = 1, bool justplot = 1)
{

    if (justplot)
    {
        plotting_helper *p_1 = new plotting_helper();
        p_1->setTDRStyle();
        TFile *f1 = new TFile("./rootfile/Fake_PbPb_new.root", "READ");

        TH1D *muon_pT = (TH1D *)f1->Get("muon_pT");
        TH1D *muon_eta = (TH1D *)f1->Get("muon_eta");
        TH1D *muon_phi = (TH1D *)f1->Get("muon_phi");

        TH1D *Fake_FA_nominal_1_1 = (TH1D *)f1->Get("Fake_FA_nominal_1_1");
        TH1D *Fake_FA_nominal_2_2 = (TH1D *)f1->Get("Fake_FA_nominal_2_2");
        TH1D *Fake_FA_nominal_22_22 = (TH1D *)f1->Get("Fake_FA_nominal_22_22");

        TH1D *Fake_FA_nominal_1_0 = (TH1D *)f1->Get("Fake_FA_nominal_1_0");
        TH1D *Fake_FA_nominal_2_1_0 = (TH1D *)f1->Get("Fake_FA_nominal_2_1_0");
        TH1D *Fake_FA_nominal_22_1_0 = (TH1D *)f1->Get("Fake_FA_nominal_22_1_0");

        TH1D *Fake_FA_nominal_1_1_pT = (TH1D *)f1->Get("Fake_FA_nominal_1_1_pT");
        TH1D *Fake_FA_nominal_2_2_pT = (TH1D *)f1->Get("Fake_FA_nominal_2_2_pT");
        TH1D *Fake_FA_nominal_22_22_pT = (TH1D *)f1->Get("Fake_FA_nominal_22_22_pT");

        TH1D *Fake_FA_nominal_1_0_pT = (TH1D *)f1->Get("Fake_FA_nominal_1_0_pT");
        TH1D *Fake_FA_nominal_1_0_pT_all_muon = (TH1D *)f1->Get("Fake_FA_nominal_1_0_pT_all_muon");
        TH1D *Fake_FA_nominal_2_1_0_pT = (TH1D *)f1->Get("Fake_FA_nominal_2_1_0_pT");
        TH1D *Fake_FA_nominal_22_1_0_pT = (TH1D *)f1->Get("Fake_FA_nominal_22_1_0_pT");

        TH1D *Fake_FA_nominal_1_1_eta = (TH1D *)f1->Get("Fake_FA_nominal_1_1_eta");
        TH1D *Fake_FA_nominal_2_2_eta = (TH1D *)f1->Get("Fake_FA_nominal_2_2_eta");
        TH1D *Fake_FA_nominal_22_22_eta = (TH1D *)f1->Get("Fake_FA_nominal_22_22_eta");

        TH1D *Fake_FA_nominal_1_0_eta = (TH1D *)f1->Get("Fake_FA_nominal_1_0_eta");
        TH1D *Fake_FA_nominal_1_0_eta_all_muon = (TH1D *)f1->Get("Fake_FA_nominal_1_0_eta_all_muon");
        TH1D *Fake_FA_nominal_2_1_0_eta = (TH1D *)f1->Get("Fake_FA_nominal_2_1_0_eta");
        TH1D *Fake_FA_nominal_22_1_0_eta = (TH1D *)f1->Get("Fake_FA_nominal_22_1_0_eta");

        TH1D *Fake_FA_nominal_1_1_phi = (TH1D *)f1->Get("Fake_FA_nominal_1_1_phi");
        TH1D *Fake_FA_nominal_2_2_phi = (TH1D *)f1->Get("Fake_FA_nominal_2_2_phi");
        TH1D *Fake_FA_nominal_22_22_phi = (TH1D *)f1->Get("Fake_FA_nominal_22_22_phi");

        TH1D *Fake_FA_nominal_1_0_phi = (TH1D *)f1->Get("Fake_FA_nominal_1_0_phi");
        TH1D *Fake_FA_nominal_1_0_phi_all_muon = (TH1D *)f1->Get("Fake_FA_nominal_1_0_phi_all_muon");
        TH1D *Fake_FA_nominal_2_1_0_phi = (TH1D *)f1->Get("Fake_FA_nominal_2_1_0_phi");
        TH1D *Fake_FA_nominal_22_1_0_phi = (TH1D *)f1->Get("Fake_FA_nominal_22_1_0_phi");

        TCanvas *c1 = new TCanvas("c1", "", 800, 800);
        c1->cd();

        plotandsave(c1, muon_pT, "pT", "./fakeplots/muon_pT.png", "All reco muon pT");
        plotandsave(c1, muon_eta, "eta", "./fakeplots/muon_eta.png", "All reco muon eta");
        plotandsave(c1, muon_phi, "phi", "./fakeplots/muon_phi.png", "All reco muon phi");

        plotandsave(c1, Fake_FA_nominal_1_1, "mass", "./fakeplots/Fake_FA_nominal_1_1.png", "1 reco, 1 gen");
        plotandsave(c1, Fake_FA_nominal_2_2, "mass", "./fakeplots/Fake_FA_nominal_2_2.png", "2 reco, 2 gen");
        plotandsave(c1, Fake_FA_nominal_22_22, "mass", "./fakeplots/Fake_FA_nominal_22_22.png", ">2 reco, >2 gen, equal size");
        plotandsave(c1, Fake_FA_nominal_1_0, "mass", "./fakeplots/Fake_FA_nominal_1_0.png", "1 reco, 0 gen");
        plotandsave(c1, Fake_FA_nominal_2_1_0, "mass", "./fakeplots/Fake_FA_nominal_2_1_0.png", "2 reco, <2 gen");
        plotandsave(c1, Fake_FA_nominal_22_1_0, "mass", "./fakeplots/Fake_FA_nominal_22_1_0.png", ">2 reco, gen < reco");

        plotandsave(c1, Fake_FA_nominal_1_1_pT, "pT", "./fakeplots/Fake_FA_nominal_1_1_pT.png", "1 reco, 1 gen");
        plotandsave(c1, Fake_FA_nominal_2_2_pT, "pT", "./fakeplots/Fake_FA_nominal_2_2_pT.png", "2 reco, 2 gen");
        plotandsave(c1, Fake_FA_nominal_22_22_pT, "pT", "./fakeplots/Fake_FA_nominal_22_22_pT.png", ">2 reco, >2 gen, equal size");
        plotandsave(c1, Fake_FA_nominal_1_0_pT, "pT", "./fakeplots/Fake_FA_nominal_1_0_pT.png", "1 reco, 0 gen");
        plotandsave(c1, Fake_FA_nominal_1_0_pT_all_muon, "pT", "./fakeplots/Fake_FA_nominal_1_0_pT_allmuon.png", "1 reco, 0 gen, all muon");
        plotandsave(c1, Fake_FA_nominal_2_1_0_pT, "pT", "./fakeplots/Fake_FA_nominal_2_1_0_pT.png", "2 reco, <2 gen");
        plotandsave(c1, Fake_FA_nominal_22_1_0_pT, "pT", "./fakeplots/Fake_FA_nominal_22_1_0_pT.png", ">2 reco, gen < reco");

        plotandsave(c1, Fake_FA_nominal_1_1_eta, "eta", "./fakeplots/Fake_FA_nominal_1_1_eta.png", "1 reco, 1 gen");
        plotandsave(c1, Fake_FA_nominal_2_2_eta, "eta", "./fakeplots/Fake_FA_nominal_2_2_eta.png", "2 reco, 2 gen");
        plotandsave(c1, Fake_FA_nominal_22_22_eta, "eta", "./fakeplots/Fake_FA_nominal_22_22_eta.png", ">2 reco, >2 gen, equal size");
        plotandsave(c1, Fake_FA_nominal_1_0_eta, "eta", "./fakeplots/Fake_FA_nominal_1_0_eta.png", "1 reco, 0 gen");
        plotandsave(c1, Fake_FA_nominal_1_0_eta_all_muon, "eta", "./fakeplots/Fake_FA_nominal_1_0_eta_allmuon.png", "1 reco, 0 gen, all muon");
        plotandsave(c1, Fake_FA_nominal_2_1_0_eta, "eta", "./fakeplots/Fake_FA_nominal_2_1_0_eta.png", "2 reco, <2 gen");
        plotandsave(c1, Fake_FA_nominal_22_1_0_eta, "eta", "./fakeplots/Fake_FA_nominal_22_1_0_eta.png", ">2 reco, gen < reco");

        plotandsave(c1, Fake_FA_nominal_1_1_phi, "phi", "./fakeplots/Fake_FA_nominal_1_1_phi.png", "1 reco, 1 gen");
        plotandsave(c1, Fake_FA_nominal_2_2_phi, "phi", "./fakeplots/Fake_FA_nominal_2_2_phi.png", "2 reco, 2 gen");
        plotandsave(c1, Fake_FA_nominal_22_22_phi, "phi", "./fakeplots/Fake_FA_nominal_22_22_phi.png", ">2 reco, >2 gen, equal size");
        plotandsave(c1, Fake_FA_nominal_1_0_phi, "phi", "./fakeplots/Fake_FA_nominal_1_0_phi.png", "1 reco, 0 gen");
        plotandsave(c1, Fake_FA_nominal_1_0_phi_all_muon, "phi", "./fakeplots/Fake_FA_nominal_1_0_phi_allmuon.png", "1 reco, 0 gen, all muon");
        plotandsave(c1, Fake_FA_nominal_2_1_0_phi, "phi", "./fakeplots/Fake_FA_nominal_2_1_0_phi.png", "2 reco, <2 gen");
        plotandsave(c1, Fake_FA_nominal_22_1_0_phi, "phi", "./fakeplots/Fake_FA_nominal_22_1_0_phi.png", ">2 reco, gen < reco");
    }
    else
    {

        MC_18 *s = new MC_18();

        auto start = std::chrono::high_resolution_clock::now();
        TH1::SetDefaultSumw2();

        const int nbins_mass_shift = 42;
        const int nbins_smear = 42;
        const int nbins_cent = 11;

        TH1D *pt_gen = new TH1D("pt_gen", "", 100, 0, 200);
        pt_gen->SetDrawOption("HIST");

        TH1D *muon_pT = new TH1D("muon_pT", "", 100, 0, 200);
        TH1D *muon_eta = new TH1D("muon_eta", "", 20, -2.4, 2.4);
        TH1D *muon_phi = new TH1D("muon_phi", "", 20, -3.14, 3.14);

        TH1D *Fake_FA_nominal_1_1 = new TH1D("Fake_FA_nominal_1_1", "", 60, 60, 120);
        TH1D *Fake_FA_nominal_2_2 = new TH1D("Fake_FA_nominal_2_2", "", 60, 60, 120);
        TH1D *Fake_FA_nominal_22_22 = new TH1D("Fake_FA_nominal_22_22", "", 60, 60, 120);

        TH1D *Fake_FA_nominal_1_0 = new TH1D("Fake_FA_nominal_1_0", "", 60, 60, 120);
        TH1D *Fake_FA_nominal_2_1_0 = new TH1D("Fake_FA_nominal_2_1_0", "", 60, 60, 120);
        TH1D *Fake_FA_nominal_22_1_0 = new TH1D("Fake_FA_nominal_22_1_0", "", 60, 60, 120);

        TH1D *Fake_FA_nominal_1_1_pT = new TH1D("Fake_FA_nominal_1_1_pT", "", 50, 0, 200);
        TH1D *Fake_FA_nominal_2_2_pT = new TH1D("Fake_FA_nominal_2_2_pT", "", 50, 0, 200);
        TH1D *Fake_FA_nominal_22_22_pT = new TH1D("Fake_FA_nominal_22_22_pT", "", 50, 0, 200);

        TH1D *Fake_FA_nominal_1_0_pT = new TH1D("Fake_FA_nominal_1_0_pT", "", 50, 0, 200);
        TH1D *Fake_FA_nominal_1_0_pT_all_muon = new TH1D("Fake_FA_nominal_1_0_pT_all_muon", "", 50, 0, 200);
        TH1D *Fake_FA_nominal_2_1_0_pT = new TH1D("Fake_FA_nominal_2_1_0_pT", "", 50, 0, 200);
        TH1D *Fake_FA_nominal_22_1_0_pT = new TH1D("Fake_FA_nominal_22_1_0_pT", "", 50, 0, 200);

        TH1D *Fake_FA_nominal_1_1_eta = new TH1D("Fake_FA_nominal_1_1_eta", "", 20, -2.4, 2.4);
        TH1D *Fake_FA_nominal_2_2_eta = new TH1D("Fake_FA_nominal_2_2_eta", "", 20, -2.4, 2.4);
        TH1D *Fake_FA_nominal_22_22_eta = new TH1D("Fake_FA_nominal_22_22_eta", "", 20, -2.4, 2.4);

        TH1D *Fake_FA_nominal_1_0_eta = new TH1D("Fake_FA_nominal_1_0_eta", "", 20, -2.4, 2.4);
        TH1D *Fake_FA_nominal_1_0_eta_all_muon = new TH1D("Fake_FA_nominal_1_0_eta_all_muon", "", 20, -2.4, 2.4);
        TH1D *Fake_FA_nominal_2_1_0_eta = new TH1D("Fake_FA_nominal_2_1_0_eta", "", 20, -2.4, 2.4);
        TH1D *Fake_FA_nominal_22_1_0_eta = new TH1D("Fake_FA_nominal_22_1_0_eta", "", 20, -2.4, 2.4);

        TH1D *Fake_FA_nominal_1_1_phi = new TH1D("Fake_FA_nominal_1_1_phi", "", 20, -3.14, 3.14);
        TH1D *Fake_FA_nominal_2_2_phi = new TH1D("Fake_FA_nominal_2_2_phi", "", 20, -3.14, 3.14);
        TH1D *Fake_FA_nominal_22_22_phi = new TH1D("Fake_FA_nominal_22_22_phi", "", 20, -3.14, 3.14);

        TH1D *Fake_FA_nominal_1_0_phi = new TH1D("Fake_FA_nominal_1_0_phi", "", 20, -3.14, 3.14);
        TH1D *Fake_FA_nominal_1_0_phi_all_muon = new TH1D("Fake_FA_nominal_1_0_phi_all_muon", "", 20, -3.14, 3.14);
        TH1D *Fake_FA_nominal_2_1_0_phi = new TH1D("Fake_FA_nominal_2_1_0_phi", "", 20, -3.14, 3.14);
        TH1D *Fake_FA_nominal_22_1_0_phi = new TH1D("Fake_FA_nominal_22_1_0_phi", "", 20, -3.14, 3.14);

        gSystem->Load("./header/libDict.so");

        MCReweight *vzRW = new MCReweight("WeightsAndEfficiencies_Z2mummu/vzReweight.root", "WeightsAndEfficiencies_Z2mummu/centralityFlatteningWeight.root");
        PtReweightSpectrum *spectrumRW = new PtReweightSpectrum("WeightsAndEfficiencies_Z2mummu/ptSpectrumReweighting.root");
        MuonTnP *tnp = new MuonTnP();

        s->SetupRootfile(1, 1);
        s->SetupBranches(0);

        TEfficiency *e[11];
        TEfficiency *e_acooff[11];

        TFile *eff_f1 = new TFile("./rootfile/mc_eff.root", "READ");
        for (int i = 0; i < 11; i++)
        {
            e[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
            e_acooff[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_withAco_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
        }

        // Histogram setup
        TH1D *event_weight = new TH1D("event_weight", "event_weight", 1, 0, 1);

        cout << "entires is " << s->t1->GetEntries() << endl;

        Long64_t totalEntries = s->t1->GetEntries();

        int counter1_1 = 0;
        int counter2_2 = 0;
        int counter22_22 = 0;

        int counter1_0 = 0;
        int counter2_1_0 = 0;
        int counter22_1_0 = 0;

        int counter1_1_counterduplicate = 0;
        int counter2_2_counterduplicate = 0;
        int counter22_22_counterduplicate = 0;

        int counter1_0_counterduplicate = 0;
        int counter2_1_0_counterduplicate = 0;
        int counter22_1_0_counterduplicate = 0;

        int counter_gen_more_than_reco = 0;

        int eventcounter = 0;

        bool testbool = true;

        for (int i = 0; i < totalEntries; i++)
        {
            double percentage = 100.0 * (i - 0) / (totalEntries - 0);

            if ((i - 0) % 1000 == 0)
            {

                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = now - start;

                // Estimate remaining time based on entries processed within the specified range
                double remaining_time = elapsed.count() * (totalEntries - i) / (i - 0 + 1);

                // Convert remaining time to a human-readable format (hours, minutes, seconds)
                int hours = static_cast<int>(remaining_time) / 3600;
                int minutes = (static_cast<int>(remaining_time) % 3600) / 60;
                int seconds = static_cast<int>(remaining_time) % 60;

                // Output progress and estimated time remaining
                std::cout << "Progress: " << percentage << "% completed, "
                          << "Estimated time remaining: "
                          << hours << "h " << minutes << "m " << seconds << "s\r"
                          << std::flush;
            }

            // W part
            s->t1->GetEntry(i);

            // event selection cut
            // hfCoincFilter2Th4 = 1,
            // primaryVertexFilter = 2,
            // clusterCompatibilityFilter = 3,

            if (!(s->evtSel[1]))
                continue;
            if (!(s->evtSel[2]))
                continue;
            if (!(s->evtSel[3]))
                continue;
            if (abs(s->bestvtxZ) > 15)
                continue;

            int hiBin = 0;
            hiBin = s->getcentrality(3);

            float eventweight = 1.0;
            eventweight = vzRW->reweightFactor(s->bestvtxZ) * vzRW->reweightFactorCent(hiBin) * s->Ncoll[hiBin] * (s->weightLHE_gen->at(1080) / 10000.0);
            // cout << eventweight << endl;
            event_weight->Fill(0.5, eventweight);

            if (!(s->trigHLT[6]))
                continue;

            bool isTau = 0;
            for (unsigned int j = 0; j < s->candSize_gen; j++)
            {
                if (TMath::Abs(s->DecayID_gen[j]) == 23)
                {
                    double ptWeight = spectrumRW->getReweightFactorMuon(s->pT_gen[j]);
                    // cout << "ptWeight is " << ptWeight << endl;
                    eventweight *= ptWeight;
                    isTau = false;
                    break;
                }
                if (TMath::Abs(s->DecayID_gen[j]) == 15)
                {
                    isTau = 1;
                    break;
                }
            }

            for (int gen = 0; gen < s->candSize_gen; gen++)
            {
                double gen_1 = s->pTD1_gen[gen];
                double gen_2 = s->pTD2_gen[gen];

                if (gen_1 < gen_2)
                {
                    pt_gen->Fill(gen_1, eventweight);
                }
                else
                {
                    pt_gen->Fill(gen_2, eventweight);
                }
            }

            // This is for reco level
            for (unsigned int j = 0; j < s->candSize; j++)
            {
                if (s->mass[j] < s->masslowlimit || s->mass[j] > s->masshighlimit)
                    continue;
                if (s->pTD1[j] < s->ptlowlimit)
                    continue;
                if (s->pTD2[j] < s->ptlowlimit)
                    continue;
                if (abs(s->EtaD1[j]) > 2.4)
                    continue;
                if (abs(s->EtaD2[j]) > 2.4)
                    continue;
                if (s->VtxProb[j] < 0.001)
                    continue;
                if (abs(s->y[j]) > 2.4)
                    continue;
                if (!(s->tightMuon1Cut(j, "POG") && s->tightMuon2Cut(j, "POG")))
                    continue;

                bool isDaughter1Trigger = s->trigMuon1->at(6).at(j); // 6 is HLT_HIL3Mu12;
                bool isDaughter2Trigger = s->trigMuon2->at(6).at(j);
                if (!(isDaughter1Trigger || isDaughter2Trigger))
                    continue;

                bool isOppositeSign = s->chargeD1[j] != s->chargeD2[j];

                if (isOppositeSign)
                {
                    bool isgenmatching = false;
                    double matchedgenindex = -99;
                    for (int gen_index = 0; gen_index < s->candSize_gen; gen_index++)
                    {
                        if (s->RecIdx_gen[gen_index] == j)
                        {
                            isgenmatching = true;
                            matchedgenindex = gen_index;
                            break;
                        }
                    }

                    if (!isTau)
                    {
                        double efficiency = s->getEfficiency(e[10], s->y[j], s->pT[j]);
                        muon_pT->Fill(s->pTD1[j], 1.0 / efficiency * eventweight);
                        muon_eta->Fill(s->EtaD1[j], 1.0 / efficiency * eventweight);
                        muon_phi->Fill(s->PhiD1[j], 1.0 / efficiency * eventweight);

                        muon_pT->Fill(s->pTD2[j], 1.0 / efficiency * eventweight);
                        muon_eta->Fill(s->EtaD2[j], 1.0 / efficiency * eventweight);
                        muon_phi->Fill(s->PhiD2[j], 1.0 / efficiency * eventweight);
                        // this is for certain reco candidate, I know if it has gen-reco matching, I know if it pass kinematic cuts, I can compare it to other reco muon as well.
                        if (!isgenmatching)
                        {
                            if (j == 0)
                            {
                                // The counter here only works for event level, so ++ only for the 0 candidate in the event
                                eventcounter++;
                            }

                            if (s->candSize == s->candSize_gen)
                            {
                                if (s->candSize == 1)
                                {
                                    counter1_1++;
                                    Fake_FA_nominal_1_1->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                    if (checkduplicate(j, s))
                                    {
                                        counter1_1_counterduplicate++;
                                        Fake_FA_nominal_1_1_pT->Fill(outputduplicate(j, s, "pT"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_1_1_eta->Fill(outputduplicate(j, s, "eta"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_1_1_phi->Fill(outputduplicate(j, s, "phi"), 1.0 / efficiency * eventweight);
                                    }
                                }
                                else if (s->candSize == 2)
                                {
                                    counter2_2++;
                                    Fake_FA_nominal_2_2->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                    if (checkduplicate(j, s))
                                    {
                                        counter2_2_counterduplicate++;
                                        Fake_FA_nominal_2_2_pT->Fill(outputduplicate(j, s, "pT"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_2_2_eta->Fill(outputduplicate(j, s, "eta"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_2_2_phi->Fill(outputduplicate(j, s, "phi"), 1.0 / efficiency * eventweight);
                                    }
                                }
                                else
                                {
                                    counter22_22++;
                                    Fake_FA_nominal_22_22->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                    if (checkduplicate(j, s))
                                    {
                                        counter22_22_counterduplicate++;
                                        Fake_FA_nominal_22_22_pT->Fill(outputduplicate(j, s, "pT"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_22_22_eta->Fill(outputduplicate(j, s, "eta"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_22_22_phi->Fill(outputduplicate(j, s, "phi"), 1.0 / efficiency * eventweight);
                                    }
                                }
                            }
                            else
                            {
                                if (s->candSize < s->candSize_gen)
                                {
                                    if (j == 0)
                                    {
                                        counter_gen_more_than_reco++;
                                        cout << "You should never see this, this means candSize < gen size " << endl;
                                    }
                                }

                                if (s->candSize == 1)
                                {
                                    counter1_0++;
                                    Fake_FA_nominal_1_0->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_nominal_1_0_pT_all_muon->Fill(s->pTD1[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_nominal_1_0_pT_all_muon->Fill(s->pTD2[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_nominal_1_0_eta_all_muon->Fill(s->EtaD1[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_nominal_1_0_eta_all_muon->Fill(s->EtaD2[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_nominal_1_0_phi_all_muon->Fill(s->PhiD1[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_nominal_1_0_phi_all_muon->Fill(s->PhiD2[j], 1.0 / efficiency * eventweight);
                                    if (checkduplicate(j, s))
                                    {
                                        counter1_0_counterduplicate++;
                                        Fake_FA_nominal_1_0_pT->Fill(outputduplicate(j, s, "pT"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_1_0_eta->Fill(outputduplicate(j, s, "eta"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_1_0_phi->Fill(outputduplicate(j, s, "phi"), 1.0 / efficiency * eventweight);
                                    }
                                }
                                else if (s->candSize == 2)
                                {
                                    counter2_1_0++;
                                    Fake_FA_nominal_2_1_0->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                    if (checkduplicate(j, s))
                                    {
                                        counter2_1_0_counterduplicate++;
                                        Fake_FA_nominal_2_1_0_pT->Fill(outputduplicate(j, s, "pT"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_2_1_0_eta->Fill(outputduplicate(j, s, "eta"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_2_1_0_phi->Fill(outputduplicate(j, s, "phi"), 1.0 / efficiency * eventweight);
                                    }
                                }
                                else
                                {
                                    counter22_1_0++;
                                    Fake_FA_nominal_22_1_0->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                    if (checkduplicate(j, s))
                                    {
                                        counter22_1_0_counterduplicate++;
                                        Fake_FA_nominal_22_1_0_pT->Fill(outputduplicate(j, s, "pT"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_22_1_0_eta->Fill(outputduplicate(j, s, "eta"), 1.0 / efficiency * eventweight);
                                        Fake_FA_nominal_22_1_0_phi->Fill(outputduplicate(j, s, "phi"), 1.0 / efficiency * eventweight);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        TFile *histogram_file;

        histogram_file = new TFile("./rootfile/Fake_PbPb_new.root", "UPDATE");

        histogram_file->cd();
        Fake_FA_nominal_1_1->Write("", 2);
        Fake_FA_nominal_2_2->Write("", 2);
        Fake_FA_nominal_22_22->Write("", 2);
        Fake_FA_nominal_1_0->Write("", 2);
        Fake_FA_nominal_2_1_0->Write("", 2);
        Fake_FA_nominal_22_1_0->Write("", 2);
        pt_gen->Write("", 2);

        Fake_FA_nominal_1_1_pT->Write("", 2);
        Fake_FA_nominal_2_2_pT->Write("", 2);
        Fake_FA_nominal_22_22_pT->Write("", 2);
        Fake_FA_nominal_1_0_pT->Write("", 2);
        Fake_FA_nominal_1_0_pT_all_muon->Write("", 2);
        Fake_FA_nominal_2_1_0_pT->Write("", 2);
        Fake_FA_nominal_22_1_0_pT->Write("", 2);

        Fake_FA_nominal_1_1_eta->Write("", 2);
        Fake_FA_nominal_2_2_eta->Write("", 2);
        Fake_FA_nominal_22_22_eta->Write("", 2);
        Fake_FA_nominal_1_0_eta->Write("", 2);
        Fake_FA_nominal_1_0_eta_all_muon->Write("", 2);
        Fake_FA_nominal_2_1_0_eta->Write("", 2);
        Fake_FA_nominal_22_1_0_eta->Write("", 2);

        Fake_FA_nominal_1_1_phi->Write("", 2);
        Fake_FA_nominal_2_2_phi->Write("", 2);
        Fake_FA_nominal_22_22_phi->Write("", 2);
        Fake_FA_nominal_1_0_phi->Write("", 2);
        Fake_FA_nominal_1_0_phi_all_muon->Write("", 2);
        Fake_FA_nominal_2_1_0_phi->Write("", 2);
        Fake_FA_nominal_22_1_0_phi->Write("", 2);

        muon_pT->Write("", 2);
        muon_eta->Write("", 2);
        muon_phi->Write("", 2);

        cout << "Final Stats: " << endl;
        cout << "We have " << eventcounter << " Events with fakes" << endl;
        cout << "We have " << counter_gen_more_than_reco << " gen more than reco " << endl;

        cout << "counter 1_1 = " << counter1_1 << endl;
        cout << "counter1_1_counterduplicate = " << counter1_1_counterduplicate << endl;
        cout << "counter 2_2 = " << counter2_2 << endl;
        cout << "counter2_2_counterduplicate = " << counter2_2_counterduplicate << endl;
        cout << "counter 22_22 = " << counter22_22 << endl;
        cout << "counter22_22_counterduplicate = " << counter22_22_counterduplicate << endl;

        cout << "counter 1_0 = " << counter1_0 << endl;
        cout << "counter1_0_counterduplicate = " << counter1_0_counterduplicate << endl;
        cout << "counter 2_1_0 = " << counter2_1_0 << endl;
        cout << "counter2_1_0_counterduplicate = " << counter2_1_0_counterduplicate << endl;
        cout << "counter 22_1_0 = " << counter22_1_0 << endl;
        cout << "counter22_1_0_counterduplicate = " << counter22_1_0_counterduplicate << endl;
    }
}
