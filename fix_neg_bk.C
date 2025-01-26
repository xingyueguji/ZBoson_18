void fix_neg_bk()
{
    TString Filepath = "./rootfile/normalized/";
    TString Filename[4] = {"efffile.root", "etacut_eff_file.root", "etacut_file.root", "rawfile.root"};

    for (int j = 0; j < 4; j++)
    {
        TString completename = Filepath + Filename[j];
        cout << "Accessing " << Filename[j] << endl;
        
        TFile *f1 = new TFile(completename, "UPDATE");

        for (int i = 0; i <= 10; ++i)
        {
            // Form the histogram name
            TString histName = TString::Format("Normalized_mc_bk_%d", i);

            // Retrieve the histogram
            TH1D *hist = dynamic_cast<TH1D *>(f1->Get(histName));
            if (!hist)
            {
                std::cerr << "Error: Histogram " << histName << " not found!" << std::endl;
                continue;
            }

            // Loop over the bins and set negative bin content to zero
            for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
            {
                double binContent = hist->GetBinContent(bin);
                if (binContent < 0)
                {
                    cout << "Fixing" << endl;
                    hist->SetBinContent(bin, 0);
                    hist->SetBinError(bin, 0);
                }
            }

            // Write the modified histogram back to the file
            hist->Write(histName, TObject::kOverwrite);
        }

        // Close the file
        f1->Close();
        delete f1;
    }
    std::cout << "All histograms processed and saved." << std::endl;
}