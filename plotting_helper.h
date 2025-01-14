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

class plotting_helper
{
public:
  plotting_helper();
  void areanormalize(TH1D *h_1);
  void areanormalize_TH1(TH1 *h_1);
  void luminormalize(TH1D *h_1, int opt, double weight);
  void compositeplot(TH1D *h_1, TH1D *h_2, TH1D *h_3, TH1D *h_4, TH1D *h_5, TH1D *h_6, int x, int opt);
  void acoplot(TH1D *h_1, TH1D *h_2, int x, bool rapiditycut);
  void savehistogram(TH1D *h_1, TH1D *h_2, TH1D *h_3, int x, TFile *f1);
  void setTDRStyle();

  double ttbar_XS = 69.0;
  double Wjet_XS = 21159;
  double DY_XS = 2010;

  double muLumi = 1682.8;
  double eLumi = 1670.7;
  double netLumi = 1695.6;

  double Nmb = 11775759052;

  float crossSectionModifier = 0.92623216;
  Int_t cenlowlimit[11] = {0, 10, 20, 30, 30, 0, 15, 50, 0, 14, 0};
  Int_t cenhighlimit[11] = {10, 20, 30, 100, 50, 15, 100, 100, 14, 100, 100};
};

plotting_helper::plotting_helper()
{
}

void plotting_helper::areanormalize(TH1D *h_1)
{
  double normalization_factor = h_1->Integral("width");
  h_1->Scale(1 / normalization_factor);
}

void plotting_helper::areanormalize_TH1(TH1 *h_1)
{
  double normalization_factor = h_1->Integral("width");
  h_1->Scale(1 / normalization_factor);
}

void plotting_helper::luminormalize(TH1D *h_1, int opt, double weight)
{
  // opt 1 == DY, opt 2 == W, opt 3 = TT
  if (opt == 1)
    h_1->Scale(DY_XS * crossSectionModifier / weight);
  if (opt == 2)
    h_1->Scale(Wjet_XS * crossSectionModifier / weight);
  if (opt == 3)
    h_1->Scale(ttbar_XS * crossSectionModifier / weight);
}

void plotting_helper::compositeplot(TH1D *h_1, TH1D *h_2, TH1D *h_3, TH1D *h_4, TH1D *h_5, TH1D *h_6, int x, int opt)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TCanvas *c1 = new TCanvas("", "", 600, 600);
  // c1->SetMargin(0,0,0,0);
  // TPad *pad = new TPad("pad", "pad", 0.15, 0.15, 0.9, 0.9);
  // pad->SetMargin(0.1,0.1,0.1,0.1);

  c1->SetLogy();
  c1->cd();
  // pad->Draw();
  // pad->cd();
  h_1->SetMarkerStyle(20); // Marker style (circle)
  h_1->SetMarkerSize(0.5); // Marker size

  h_1->SetMarkerColor(kRed);
  h_2->SetFillColor(kOrange + 1);
  h_3->SetFillColor(kBlue);
  h_4->SetFillColor(kGreen);
  h_5->SetFillColor(kMagenta + 2);
  h_6->SetFillColor(kGray);

  h_2->SetOption("HIST");
  h_3->SetOption("HIST");
  h_4->SetOption("HIST");
  h_5->SetOption("HIST");
  h_6->SetOption("HIST");

  for (int i = 1; i <= 30; i++)
  {
    if (h_2->GetBinContent(i) < 0)
      cout << "Negative in h_2!" << endl;
    if (h_3->GetBinContent(i) < 0)
      cout << "Negative in h_3!" << endl;
    if (h_4->GetBinContent(i) < 0)
      cout << "Negative in h_4!" << endl;
    if (h_5->GetBinContent(i) < 0)
      cout << "Negative in h_5!" << endl;
    if (h_6->GetBinContent(i) < 0)
      cout << "Negative in h_6!" << endl;
  }

  TLegend *legend = new TLegend(0.66, 0.7, 0.93, 0.93);
  legend->SetTextSize(0.03);
  legend->SetTextFont(42);
  legend->SetFillColorAlpha(kWhite, 0);
  legend->AddEntry(h_1, Form("Data-EM (%i-%i%)", cenlowlimit[x], cenhighlimit[x]), "p");
  legend->AddEntry(h_2, "Z \\rightarrow \\mu^{+} \\mu^{-}", "f");
  legend->AddEntry(h_3, "Same_sign(QCD)", "f");
  legend->AddEntry(h_4, "Z \\rightarrow \\tau^{+} \\tau^{-}", "f");
  legend->AddEntry(h_5, "W^{\\pm}", "f");
  legend->AddEntry(h_6, "t\\bar{t}", "f");
  // legend->SetBorderSize(0); // Set border size of the legend
  THStack *s1 = new THStack("s1", "");
  s1->SetTitle("");
  s1->SetMinimum(1e-5);
  s1->SetMaximum(2e-1);
  // s1->Add(h_7,"HIST");
  s1->Add(h_6, "HIST");
  s1->Add(h_5, "HIST");
  s1->Add(h_4, "HIST");
  s1->Add(h_3, "HIST");
  s1->Add(h_2, "HIST");
  s1->Draw();
  s1->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
  s1->GetXaxis()->CenterTitle();
  s1->GetYaxis()->SetTitle("Events / 4.0 GeV");
  s1->GetYaxis()->CenterTitle();
  h_1->Draw("SAME");
  legend->Draw("SAME");

  TLatex runinfo;
  runinfo.SetNDC();
  runinfo.SetTextSize(0.04);
  runinfo.SetTextFont(42);
  runinfo.DrawLatex(0.58, 0.96, "1.8 nb^{-1} (5.02 TeV PbPb)");

  TLatex logo;
  logo.SetNDC();
  logo.SetTextSize(0.04);
  logo.SetTextFont(62);
  // logo.DrawLatex(0.2,0.88,"CMS");

  TLatex logoextra;
  logoextra.SetNDC();
  logoextra.SetTextSize(0.04);
  logoextra.SetTextFont(52);
  // logoextra.DrawLatex(0.2,0.83,"Preliminary");
  if (opt == 1)
    c1->SaveAs(Form("./newcomposite/raw/composite_%i.pdf", x));
  if (opt == 2)
    c1->SaveAs(Form("./newcomposite/y/composite_%i.pdf", x));
  if (opt == 3)
    c1->SaveAs(Form("./newcomposite/eff/composite_%i.pdf", x));
  if (opt == 4)
    c1->SaveAs(Form("./newcomposite/y_eff/composite_%i.pdf", x));
  if (opt == 5)
    c1->SaveAs(Form("./newcomposite/eta/composite_%i.pdf", x));
  if (opt == 6)
    c1->SaveAs(Form("./newcomposite/eta_eff/composite_%i.pdf", x));
}
void plotting_helper::acoplot(TH1D *h_1, TH1D *h_2, int x, bool rapiditycut)
{
  TCanvas *c2 = new TCanvas("", "", 800, 600);
  c2->cd();
  h_1->SetLineColor(kBlue);
  h_2->SetLineColor(kRed);
  h_1->GetXaxis()->SetRangeUser(0, 0.1);
  h_2->GetXaxis()->SetRangeUser(0, 0.1);
  h_2->GetXaxis()->SetTitle("Acoplanarity");
  h_2->Draw("P");      // data
  h_1->Draw("P SAME"); // mc

  if (rapiditycut)
    c2->SaveAs(Form("./acoplot/A_%i_withycut.pdf", x));
  if (!rapiditycut)
    c2->SaveAs(Form("./acoplot/A_%i.pdf", x));
}

void plotting_helper::savehistogram(TH1D *h_1, TH1D *h_2, TH1D *h_3, int x, TFile *f1)
{
  f1->cd();
  h_1->Write(Form("Normalized_mc_%i", x), 2);
  h_2->Write(Form("Normalized_data_%i", x), 2);
  h_3->Write(Form("Normalized_mc_bk_%i", x), 2);
}

void plotting_helper::setTDRStyle()
{
  TStyle *tdrStyle = new TStyle("tdrStyle", "Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0); // with or without line/dashed line as border
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); // Height of canvas
  tdrStyle->SetCanvasDefW(600); // Width of canvas
  tdrStyle->SetCanvasDefX(0);   // POsition on screen (GUI related)
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo: (I want to specify it later)
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  // tdrStyle->SetHistLineColor(1);
  // tdrStyle->SetHistLineStyle(0);
  // tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  // tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  // For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g"); // 5 in width, 4 in digits after decimal, g = general
  tdrStyle->SetFuncColor(2);      // All for TF1
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  // For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  /*tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);*/
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.04);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.05, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.1);
  tdrStyle->SetTitleYOffset(1.45);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20., 20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();
}
