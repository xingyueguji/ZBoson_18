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

// C++ stuff
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

void get_all_bk_mc_fake(int type = 1)
{

    MC_18 *s = new MC_18();

    auto start = std::chrono::high_resolution_clock::now();
    TH1::SetDefaultSumw2();

    const int nbins_mass_shift = 42;
    const int nbins_smear = 42;
    const int nbins_cent = 11;

    TH1D *Fake_FA_nominal[nbins_cent];

    TH1D *Fake_FA_pT[nbins_cent];
    TH1D *Fake_FA_Eta[nbins_cent];
    TH1D *Fake_FA_Phi[nbins_cent];
    TH1D *Fake_FA_y[nbins_cent];

    TH1D *Normal_FA_pT[nbins_cent];
    TH1D *Normal_FA_Eta[nbins_cent];
    TH1D *Normal_FA_Phi[nbins_cent];
    TH1D *Normal_FA_y[nbins_cent];

    for (int i = 0; i < nbins_cent; i++)
    {
        Fake_FA_nominal[i] = new TH1D(Form("Fake_FA_nominal_%i", i), "", 120, 60, 120);

        Fake_FA_pT[i] = new TH1D(Form("Fake_FA_pT_%i", i), "", 50, 0, 200);
        Fake_FA_Eta[i] = new TH1D(Form("Fake_FA_Eta_%i", i), "", 20, -2.4, 2.4);
        Fake_FA_Phi[i] = new TH1D(Form("Fake_FA_Phi_%i", i), "", 20, -3.14, 3.14);
        Fake_FA_y[i] = new TH1D(Form("Fake_FA_y_%i", i), "", 20, -3.14, 3.14);

        Normal_FA_pT[i] = new TH1D(Form("Normal_FA_pT_%i", i), "", 50, 0, 200);
        Normal_FA_Eta[i] = new TH1D(Form("Normal_FA_Eta_%i", i), "", 20, -2.4, 2.4);
        Normal_FA_Phi[i] = new TH1D(Form("Normal_FA_Phi_%i", i), "", 20, -3.14, 3.14);
        Normal_FA_y[i] = new TH1D(Form("Normal_FA_y_%i", i), "", 20, -3.14, 3.14);
    }

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
        e_acooff[i] = (TEfficiency *)eff_f1->Get(Form("eff_noSF_noAco_%i_%i", s->cenlowlimit[i], s->cenhighlimit[i]));
    }

    // Histogram setup
    TH1D *event_weight = new TH1D("event_weight", "event_weight", 1, 0, 1);

    cout << "entires is " << s->t1->GetEntries() << endl;

    Long64_t totalEntries = s->t1->GetEntries();

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

            int centbinpositioncounter[11] = {};
            s->CentBinSearching(centbinpositioncounter, hiBin);
            // cout << "hiBin is " << hiBin << endl;

            float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(s->PhiD1[j] - s->PhiD2[j]))) / TMath::Pi();
            bool passesAco[3] = {1, 1, 1};
            if (s->pT[j] < 1.25 && acoplanarity < 0.001)
                passesAco[0] = false;

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

                for (int k = 0; k < s->centarraysize; k++)
                {
                    if (!((k < 4) || (k == 10)))
                        continue;
                    if (centbinpositioncounter[k] != 0)
                    {
                        if (!isTau)
                        {

                            if (isgenmatching)
                            {
                                double gen_mass = s->Calc_Z_gen(matchedgenindex);
                                double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);

                                Normal_FA_pT[k]->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                Normal_FA_Eta[k]->Fill(s->eta[j], 1.0 / efficiency * eventweight);
                                Normal_FA_Phi[k]->Fill(s->phi[j], 1.0 / efficiency * eventweight);
                                Normal_FA_y[k]->Fill(s->y[j], 1.0 / efficiency * eventweight);

                                if (gen_mass <= 60 || gen_mass >= 120)
                                {
                                    Fake_FA_nominal[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_pT[k]->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_Eta[k]->Fill(s->eta[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_Phi[k]->Fill(s->phi[j], 1.0 / efficiency * eventweight);
                                    Fake_FA_y[k]->Fill(s->y[j], 1.0 / efficiency * eventweight);
                                }
                            }
                            else
                            {
                                double efficiency = s->getEfficiency(e[k], s->y[j], s->pT[j]);
                                Fake_FA_nominal[k]->Fill(s->mass[j], 1.0 / efficiency * eventweight);
                                Fake_FA_pT[k]->Fill(s->pT[j], 1.0 / efficiency * eventweight);
                                Fake_FA_Eta[k]->Fill(s->eta[j], 1.0 / efficiency * eventweight);
                                Fake_FA_Phi[k]->Fill(s->phi[j], 1.0 / efficiency * eventweight);
                                Fake_FA_y[k]->Fill(s->y[j], 1.0 / efficiency * eventweight);
                            }
                        }
                    }
                }

                if (!isgenmatching)
                {
                    continue;
                }
            }
        }
    }

    TFile *histogram_file;

    histogram_file = new TFile("./rootfile/Fake_PbPb.root", "UPDATE");

    histogram_file->cd();

    for (int i = 0; i < nbins_cent; i++)
    {
        if (!((i < 4) || (i == 10)))
            continue;

        Fake_FA_nominal[i]->Write("", 2);
        Fake_FA_pT[i]->Write("", 2);
        Fake_FA_Eta[i]->Write("", 2);
        Fake_FA_Phi[i]->Write("", 2);
        Fake_FA_y[i]->Write("", 2);
        Normal_FA_pT[i]->Write("", 2);
        Normal_FA_Eta[i]->Write("", 2);
        Normal_FA_Phi[i]->Write("", 2);
        Normal_FA_y[i]->Write("", 2);
    }
}
