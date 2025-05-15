#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <map>

void count_DecayID_gen(const char *filename = "~/ROOTFILES/mc_signal_file.root", const char *tree_name = "VertexCompositeNtuple")
{
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    TDirectoryFile *d_1 = (TDirectoryFile*)file->Get("dimucontana_mc");

    // Get the TTree
    TTree *tree = (TTree *)d_1->Get(tree_name);
    if (!tree)
    {
        std::cerr << "Error: Cannot find tree " << tree_name << " in file " << filename << std::endl;
        file->Close();
        return;
    }

    // Disable all branches first to speed up the reading process
    tree->SetBranchStatus("*", 0);

    // Enable only the needed branches
    Short_t DecayID_gen[5000] = {}; // Array for decay ID
    UInt_t candSize_gen = 0;        // Number of candidates

    tree->SetBranchStatus("DecayID_gen", 1);
    tree->SetBranchStatus("candSize_gen", 1);
    tree->SetBranchAddress("DecayID_gen", DecayID_gen);
    tree->SetBranchAddress("candSize_gen", &candSize_gen);

    std::map<Short_t, int> decayID_counts; // Map to store counts of each unique ID

    // Loop over tree entries
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i)
    {
        tree->GetEntry(i);
        if (i % 100000 == 0)
        {
            cout << "This is entry " << i << endl;
        }

        for (UInt_t j = 0; j < candSize_gen; ++j)
        {
            decayID_counts[DecayID_gen[j]]++;
        }
    }

    // Print results
    std::cout << "Decay ID Counts:\n";
    for (const auto &entry : decayID_counts)
    {
        std::cout << "Decay ID: " << entry.first << " -> Count: " << entry.second << std::endl;
    }

    // Clean up
    file->Close();
}