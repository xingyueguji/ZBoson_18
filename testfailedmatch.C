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

void testfailedmatch()
{

    MC_18 *s = new MC_18();

    gSystem->Load("./header/libDict.so");

    MCReweight *vzRW = new MCReweight("WeightsAndEfficiencies_Z2mummu/vzReweight.root", "WeightsAndEfficiencies_Z2mummu/centralityFlatteningWeight.root");
    PtReweightSpectrum *spectrumRW = new PtReweightSpectrum("WeightsAndEfficiencies_Z2mummu/ptSpectrumReweighting.root");
    MuonTnP *tnp = new MuonTnP();

    s->SetupRootfile(1, 1);
    s->SetupBranches(0);

    cout << "entires is " << s->t1->GetEntries() << endl;

    Long64_t totalEntries = s->t1->GetEntries();

    int outsidecounter = 0;
    int fakecounter = 0;

    for (int i = 0; i < totalEntries; i++)
    {

        // W part
        s->t1->GetEntry(i);
        if (i % 100000 == 0)
        {
            cout << "This is event " << i << endl;
        }

        /*if (s->candSize_gen == 0)
        {
            if (s->matchGEN[0] == true)
            {
                if (s->candSize == 0) continue;
                cout << "We have matchGEN false case " << endl;
                cout << "This is event " << i << endl;
            }
        }*/

        // Don't use matchGEN, event 39, we have candsize_gen = 0, candsize = 2 and matchGEN = 1 and 0.

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

        if (!(s->trigHLT[6]))
            continue;

        bool isTau = 0;
        for (unsigned int j = 0; j < s->candSize_gen; j++)
        {
            if (TMath::Abs(s->DecayID_gen[j]) == 23)
            {
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
                        // cout << "matchGEN result is " << s->matchGEN[j] << endl;
                        break;
                    }
                }

                if (!isTau)
                {
                    if (passesAco[0])
                    {
                        if (isgenmatching)
                        {
                            if (s->matchGEN[j] == false)
                            {
                                cout << "We have gen reco mataching disagree with match GEN " << endl;
                            }
                            double gen_mass = s->Calc_Z_gen(matchedgenindex);
                            if (gen_mass <= 60 || gen_mass >= 120)
                            {
                                outsidecounter++;
                            }
                        }
                        else
                        {
                            fakecounter++;
                            /*cout << "We Have one fake candidate" << endl;
                            cout << "The Reco cand size is " << s->candSize << endl;
                            cout << "This is the " << j << " Reco Z with matchGEN " << s->matchGEN[j] << endl;
                            cout << "The Gen cand size is " << s->candSize_gen << endl;
                            cout << "Here's all gen level info: " << endl;

                            for (int gen_index = 0; gen_index < s->candSize_gen; gen_index++)
                            {
                                cout << "The " << gen_index << " Gen particle has: " << endl;
                                cout << "RecIdx_gen = " << s->RecIdx_gen[gen_index] << " DecayID_gen = " << s->DecayID_gen[gen_index] << " MotherID_gen = " << s->MotherID_gen[gen_index] << " PID_gen = " << s->PID_gen[gen_index] << " status_gen = " << s->status_gen[gen_index] << endl;
                            }*/

                            //We have outsidecounter = 683 fakecounter = 34345
                        }
                    }
                }
            }
        }
    }
    cout << "We have outsidecounter = " << outsidecounter << " fakecounter = " << fakecounter << endl;
}