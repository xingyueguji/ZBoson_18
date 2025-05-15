#pragma once

#include "TMath.h"
// Header file for ROOT classes
#include <TROOT.h>
#include <TChain.h>
#include <TInterpreter.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TFile.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TPad.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include "TEfficiency.h"

// Header file for c++ classes
#include <iostream>
#include <string>
#include <vector>
#include <map>


typedef std::vector<std::vector<UChar_t>> UCharVecVec;

class MC_18
{
public:
	MC_18();
	~MC_18();

	void SetupRootfile(Int_t x, Int_t y);
	void SetupBranches(Int_t x);
	// bool SetupCut();
	void Plot_and_Pull(TH1D *hist, TH1D *hist2, TCanvas *c1, TPad *p1, TPad *p2);
	void Plot_and_Pull_gen(TH1D *hist, TH1D *hist2, TCanvas *c1, TPad *p1, TPad *p2);
	RooFormulaVar *CreateRatio(const RooAbsPdf *pdf1, const RooAbsPdf *pdf2, RooRealVar *x);
	Float_t Calc_Z_gen(int indx);
	void CentBinSearching(int *arr, int hiBin);
	void getrapidity(int x);
	int getcentrality(int opt);
	void forceConsistency(TH2D *top, TH2D *bot);
	Bool_t tightMuon1Cut(int x, string type);
	Bool_t tightMuon2Cut(int x, string type);
	double getEfficiency(TEfficiency *e, double y, double pt);

	// Here's lxplus version
	// TString MC_File_Path = "root://cms-xrd-global.cern.ch///store/group/phys_heavyions/abaty/Zmumu_Samples_MC/VertexCompositeTree_DYJetsToLL_TuneCP5_HydjetDrumMB_Pythia8_HINPbPbAutumn18_DiMuMC_20190808.root";
	// TString MC_W_File_Path = "root://cms-xrd-global.cern.ch///store/group/phys_heavyions/abaty/Zmumu_Samples_MC/VertexCompositeTree_WJetsToLNu_TuneCP5_HydjetDrumMB_Pythia8_HINPbPbAutumn18_DiMuMC.root";
	// TString MC_tt_File_Path = "root://eoshome-a.cern.ch//eos/user/a/abaty/Z_Datasets_vtxTrees_mumu/VertexCompositeTree_TTJets_TuneCP5_HydjetDrumMB_Pythia8_HINPbPbAutumn18_DiMuMC.root";
	// TString Data_File_Path = "root://eoshome-a.cern.ch//eos/user/a/abaty/Z_Datasets_vtxTrees_mumu/VertexCompositeTree_OppositeSignSkim_HLTbit6_gt0Cands_v2.root";
	// TString Data_File_same_Path = "root://eoshome-a.cern.ch//eos/user/a/abaty/Z_Datasets_vtxTrees_mumu/VertexCompositeTree_SameSignSkim_HLTbit6_gt0Cands_v2.root";

	// Here's local version
	TString MC_File_Path = "~/ROOTFILES/mc_signal_file.root";
	TString MC_W_File_Path = "~/ROOTFILES/VertexCompositeTree_WJetsToLNu_TuneCP5_HydjetDrumMB_Pythia8_HINPbPbAutumn18_DiMuMC.root";
	TString MC_tt_File_Path = "~/ROOTFILES/VertexCompositeTree_TTJets_TuneCP5_HydjetDrumMB_Pythia8_HINPbPbAutumn18_DiMuMC.root";
	TString Data_File_Path = "~/ROOTFILES/VertexCompositeTree_OppositeSignSkim_HLTbit6_gt0Cands_v2.root";
	TString Data_File_same_Path = "~/ROOTFILES/VertexCompositeTree_SameSignSkim_HLTbit6_gt0Cands_v2.root";

	TFile *f1;
	TDirectoryFile *dir1;
	TTree *t1;

	TFile *f2;
	TDirectoryFile *dir2;
	TTree *t2;

	// These are reco level;
	UInt_t candSize = 0;
	Short_t centrality = 0;
	Float_t HFsumETPlus = 0;
	Float_t HFsumETMinus = 0;
	Float_t bestvtxZ = 0;
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
	Int_t PIDD2[5000] = {};
	Int_t PIDD1[5000] = {};
	Bool_t evtSel[4] = {};
	Float_t y[5000] = {};
	Bool_t tightMuon1[5000] = {};
	Bool_t tightMuon2[5000] = {};

	Bool_t GlbMuon1[5000] = {};
	Bool_t GlbMuon2[5000] = {};
	Bool_t PFMuon1[5000] = {};
	Bool_t PFMuon2[5000] = {};
	Float_t GlbTrkChiD1[5000] = {};
	Float_t GlbTrkChiD2[5000] = {};
	Short_t nMuonHitD1[5000] = {};
	Short_t nMuonHitD2[5000] = {};
	Short_t nMatchedStationD1[5000] = {};
	Short_t nMatchedStationD2[5000] = {};
	Short_t nPixelHitD1[5000] = {};
	Short_t nPixelHitD2[5000] = {};
	Short_t nTrackerLayerD1[5000] = {};
	Short_t nTrackerLayerD2[5000] = {};
	Float_t muondXYD1[5000] = {};
	Float_t muondXYD2[5000] = {};
	Float_t muondZD1[5000] = {};
	Float_t muondZD2[5000] = {};
	Bool_t trigHLT[9] = {};

	vector<float> *weightLHE_gen = 0;

	UCharVecVec *trigMuon1 = 0;
	UCharVecVec *trigMuon2 = 0;

	// These are I believe gen level
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
	Int_t PID_gen[5000] = {};
	Float_t y_gen[5000] = {};
	Float_t pT_gen[5000] = {};
	Short_t RecIdx_gen[5000] = {};
	Short_t chargeD1_gen[5000] = {};
	Short_t chargeD2_gen[5000] = {};
	Int_t MotherID_gen[5000] = {};
	Short_t status_gen[5000] = {};
	Bool_t matchGEN[5000] = {};

	Int_t nentries = 0;

	// Cut information
	Int_t masslowlimit = 60;
	Int_t masshighlimit = 120;
	Double_t ptlowlimit = 20; // check 15
	Double_t etalimit = 2.4;
	Double_t vtxproblowlimit = 0.001;

	// centrality range def
	Int_t centarraysize = 11;
	Int_t cenlowlimit[11] = {0, 10, 20, 30, 30, 0, 15, 50, 0, 14, 0};
	Int_t cenhighlimit[11] = {10, 20, 30, 100, 50, 15, 100, 100, 14, 100, 100};

	Float_t rapidity = 0;

	const Double_t binTableMC_Drum5F[200 + 1] = {0, 12.2187, 13.0371, 13.7674, 14.5129, 15.2603, 16.0086, 16.7623, 17.5335, 18.3283, 19.1596, 19.9989, 20.8532, 21.7297, 22.6773, 23.6313, 24.6208, 25.6155, 26.6585, 27.7223, 28.8632, 30.041, 31.2865, 32.5431, 33.8655, 35.2539, 36.6912, 38.2064, 39.7876, 41.4818, 43.2416, 45.0605, 46.9652, 48.9918, 51.1, 53.2417, 55.5094, 57.9209, 60.3817, 62.9778, 65.6099, 68.4352, 71.3543, 74.4154, 77.6252, 80.8425, 84.1611, 87.7395, 91.3973, 95.1286, 99.0571, 103.185, 107.482, 111.929, 116.45, 121.178, 126.081, 130.995, 136.171, 141.612, 147.298, 153.139, 159.419, 165.633, 172.114, 178.881, 185.844, 192.845, 200.244, 207.83, 215.529, 223.489, 231.878, 240.254, 249.319, 258.303, 267.508, 277.037, 286.729, 296.845, 307.458, 317.882, 328.787, 340.074, 351.295, 362.979, 375.125, 387.197, 399.604, 412.516, 425.683, 439.001, 452.667, 466.816, 481.007, 495.679, 510.588, 526.138, 541.782, 557.641, 574.141, 591.071, 608.379, 626.068, 643.616, 661.885, 680.288, 699.449, 718.925, 738.968, 758.983, 779.459, 800.376, 821.638, 843.555, 865.771, 888.339, 911.031, 934.979, 958.56, 982.582, 1007.02, 1031.9, 1057.81, 1084.01, 1111.71, 1138.21, 1165.72, 1193.73, 1221.65, 1251.51, 1281.23, 1311.01, 1341.1, 1372.4, 1404.29, 1436.52, 1468.65, 1501.91, 1535.56, 1569.69, 1604.69, 1640.65, 1676.05, 1712.62, 1749.28, 1787.43, 1825.89, 1866.07, 1906.58, 1947.84, 1989.66, 2031.4, 2072.8, 2115.32, 2159.5, 2205.23, 2252.68, 2298.58, 2345.65, 2393.36, 2442.87, 2491.45, 2541.04, 2592.81, 2645.52, 2699.1, 2753.29, 2807.93, 2864.37, 2922.6, 2979.42, 3038.68, 3098.72, 3159.29, 3221.66, 3285.9, 3350.95, 3415.81, 3482.69, 3552.62, 3623.61, 3694.63, 3767.25, 3840.28, 3917.04, 3993.66, 4073.36, 4154.33, 4238.13, 4322.21, 4409.83, 4498.89, 4589.72, 4681.56, 4777.09, 4877.95, 4987.05, 5113.04, 5279.58, 6242.82};
	const float Ncoll[200] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
	const Double_t binTableData[200 + 1] = {0, 10.5072, 11.2099, 11.8364, 12.478, 13.1194, 13.7623, 14.4081, 15.0709, 15.7532, 16.4673, 17.1881, 17.923, 18.673, 19.4865, 20.3033, 21.1536, 22.0086, 22.9046, 23.8196, 24.7924, 25.8082, 26.8714, 27.9481, 29.0828, 30.2757, 31.5043, 32.8044, 34.1572, 35.6142, 37.1211, 38.6798, 40.3116, 42.0398, 43.8572, 45.6977, 47.6312, 49.6899, 51.815, 54.028, 56.3037, 58.7091, 61.2024, 63.8353, 66.5926, 69.3617, 72.2068, 75.2459, 78.3873, 81.5916, 84.9419, 88.498, 92.1789, 95.9582, 99.8431, 103.739, 107.78, 111.97, 116.312, 120.806, 125.46, 130.269, 135.247, 140.389, 145.713, 151.212, 156.871, 162.729, 168.762, 174.998, 181.424, 188.063, 194.907, 201.942, 209.19, 216.683, 224.37, 232.291, 240.43, 248.807, 257.416, 266.256, 275.348, 284.668, 294.216, 304.053, 314.142, 324.488, 335.101, 345.974, 357.116, 368.547, 380.283, 392.29, 404.564, 417.122, 429.968, 443.116, 456.577, 470.357, 484.422, 498.78, 513.473, 528.479, 543.813, 559.445, 575.411, 591.724, 608.352, 625.344, 642.686, 660.361, 678.371, 696.749, 715.485, 734.608, 754.068, 773.846, 794.046, 814.649, 835.608, 856.972, 878.719, 900.887, 923.409, 946.374, 969.674, 993.435, 1017.62, 1042.21, 1067.28, 1092.72, 1118.64, 1144.96, 1171.71, 1198.98, 1226.67, 1254.82, 1283.46, 1312.65, 1342.21, 1372.27, 1402.85, 1433.93, 1465.49, 1497.62, 1530.29, 1563.49, 1597.22, 1631.49, 1666.37, 1701.8, 1737.75, 1774.35, 1811.51, 1849.29, 1887.75, 1926.79, 1966.6, 2006.97, 2047.99, 2089.71, 2132.1, 2175.23, 2219.17, 2263.72, 2309.2, 2355.43, 2402.47, 2450.33, 2499.05, 2548.66, 2599.16, 2650.59, 2703.03, 2756.32, 2810.75, 2866.27, 2922.91, 2980.54, 3039.47, 3099.53, 3160.98, 3223.66, 3287.71, 3353.18, 3420.34, 3489.13, 3559.72, 3632.06, 3706.18, 3782.42, 3860.78, 3941.42, 4024.52, 4110.27, 4199.4, 4292.8, 4394.49, 4519.52, 5199.95};
	const Double_t binTableDataDown[200 + 1] = {0, 10.5071, 11.2094, 11.8357, 12.4763, 13.117, 13.7597, 14.4049, 15.0671, 15.7491, 16.4622, 17.1812, 17.9144, 18.6674, 19.4797, 20.2963, 21.1435, 21.9974, 22.8928, 23.8068, 24.7805, 25.7931, 26.8556, 27.9308, 29.0638, 30.2582, 31.4795, 32.7816, 34.1349, 35.5834, 37.0941, 38.6474, 40.2782, 42.0035, 43.8112, 45.6576, 47.5758, 49.6381, 51.6667, 53.7353, 55.8903, 58.1259, 60.4528, 62.8712, 65.3859, 67.9968, 70.7065, 73.5231, 76.4519, 79.4922, 82.6461, 85.9264, 89.3269, 92.8562, 96.5212, 100.322, 104.262, 108.344, 112.585, 116.971, 121.521, 126.225, 131.09, 136.127, 141.328, 146.721, 152.284, 158.014, 163.935, 170.054, 176.372, 182.878, 189.602, 196.532, 203.653, 211.017, 218.599, 226.387, 234.418, 242.667, 251.16, 259.886, 268.852, 278.071, 287.498, 297.2, 307.184, 317.409, 327.894, 338.66, 349.686, 360.996, 372.607, 384.508, 396.669, 409.133, 421.86, 434.906, 448.258, 461.916, 475.906, 490.16, 504.74, 519.663, 534.911, 550.453, 566.322, 582.525, 599.08, 615.968, 633.211, 650.805, 668.76, 687.048, 705.707, 724.774, 744.163, 763.9, 783.999, 804.528, 825.432, 846.746, 868.429, 890.523, 913.007, 935.952, 959.211, 982.919, 1007.08, 1031.63, 1056.62, 1082.08, 1107.96, 1134.24, 1160.99, 1188.22, 1215.91, 1244.06, 1272.69, 1301.85, 1331.45, 1361.51, 1392.07, 1423.18, 1454.77, 1486.93, 1519.57, 1552.81, 1586.55, 1620.87, 1655.79, 1691.26, 1727.27, 1763.93, 1801.12, 1838.97, 1877.47, 1916.61, 1956.45, 1996.89, 2038.04, 2079.84, 2122.35, 2165.52, 2209.53, 2254.24, 2299.83, 2346.19, 2393.31, 2441.28, 2490.16, 2539.86, 2590.57, 2642.16, 2694.74, 2748.23, 2802.81, 2858.47, 2915.33, 2973.2, 3032.28, 3092.56, 3154.24, 3217.19, 3281.45, 3347.18, 3414.6, 3483.65, 3554.56, 3627.2, 3701.66, 3778.25, 3856.97, 3937.98, 4021.48, 4107.62, 4197.21, 4291.05, 4393.19, 4518.6, 5199.95};
	const Double_t binTableDataUp[200 + 1] = {0, 10.5075, 11.2107, 11.838, 12.4797, 13.1213, 13.7641, 14.4124, 15.0745, 15.7577, 16.473, 17.1939, 17.9297, 18.6812, 19.4958, 20.3143, 21.1648, 22.0218, 22.9159, 23.8328, 24.8059, 25.8204, 26.89, 27.9702, 29.1042, 30.3022, 31.528, 32.8347, 34.1896, 35.6439, 37.1542, 38.7172, 40.3518, 42.091, 43.9053, 45.7415, 47.6853, 49.7457, 51.8755, 54.0983, 56.3594, 58.7848, 61.2861, 63.9228, 66.6825, 69.4421, 72.297, 75.3547, 78.4967, 81.6977, 85.0755, 88.6211, 92.3058, 96.1071, 99.9975, 104.065, 108.272, 112.512, 116.906, 121.601, 126.465, 131.482, 136.866, 142.229, 147.786, 153.546, 159.571, 165.586, 171.902, 178.419, 185.063, 191.856, 199.055, 206.261, 213.999, 221.719, 229.671, 237.84, 246.088, 254.828, 263.883, 272.907, 282.236, 291.925, 301.519, 311.477, 321.691, 332.153, 342.892, 353.878, 365.161, 376.742, 388.577, 400.684, 413.075, 425.746, 438.711, 451.989, 465.556, 479.45, 493.608, 508.077, 522.891, 538.003, 553.415, 569.151, 585.216, 601.601, 618.354, 635.422, 652.84, 670.599, 688.699, 707.161, 726.014, 745.185, 764.687, 784.557, 804.838, 825.489, 846.537, 867.951, 889.752, 911.955, 934.588, 957.52, 980.912, 1004.73, 1028.94, 1053.57, 1078.67, 1104.17, 1130.07, 1156.39, 1183.2, 1210.47, 1238.17, 1266.38, 1295.02, 1324.16, 1353.71, 1383.77, 1414.35, 1445.41, 1477, 1509.09, 1541.74, 1574.88, 1608.59, 1642.83, 1677.66, 1713.07, 1748.98, 1785.47, 1822.63, 1860.33, 1898.72, 1937.73, 1977.42, 2017.71, 2058.62, 2100.25, 2142.57, 2185.56, 2229.38, 2273.91, 2319.2, 2365.33, 2412.22, 2459.94, 2508.52, 2557.98, 2608.35, 2659.61, 2711.86, 2765, 2819.23, 2874.58, 2930.97, 2988.46, 3047.12, 3106.95, 3168.15, 3230.6, 3294.37, 3359.58, 3426.47, 3494.95, 3565.21, 3637.21, 3711.03, 3786.91, 3864.85, 3945.11, 4027.8, 4113.06, 4201.73, 4294.72, 4395.9, 4520.5, 5199.95};
};

MC_18::MC_18()
{
}

MC_18::~MC_18()
{
	f1->Close();
	delete f1;
	delete t1;
	delete dir1;
}

void MC_18::CentBinSearching(int *arr, int hiBin)
{
	for (int i = 0; i < centarraysize; i++)
	{
		if (hiBin / 2 >= cenlowlimit[i] && hiBin / 2 < cenhighlimit[i])
		{
			arr[i] = 1;
		}
		else
		{
			arr[i] = 0;
		}
	}
}

void MC_18::SetupRootfile(Int_t x, Int_t y)
{
	if (x == 1)
		f1 = TFile::Open(MC_File_Path, "READ");
	if (x == 2)
		f1 = TFile::Open(MC_W_File_Path, "READ");
	if (x == 3)
		f1 = TFile::Open(MC_tt_File_Path, "READ");
	if (x == 4)
		f1 = TFile::Open(Data_File_Path, "READ");
	if (x == 5)
		f1 = TFile::Open(Data_File_same_Path, "READ");

	if (f1->IsOpen() != 1)
	{
		if (x == 1)
			cout << "ROOTfile " << MC_File_Path << " Cannot be opened! " << endl;
		if (x == 2)
			cout << "ROOTfile " << MC_W_File_Path << " Cannot be opened! " << endl;
		if (x == 3)
			cout << "ROOTfile " << MC_tt_File_Path << " Cannot be opened! " << endl;
	}

	if (y == 1)
		dir1 = (TDirectoryFile *)f1->Get("dimucontana_mc;1");
	if (y == 2)
		dir1 = (TDirectoryFile *)f1->Get("dimucontana_wrongsign_mc;1");
	if (!dir1 && y != 0)
	{
		cout << "TDirectory dimucontana cannot be opened! " << endl;
	}
	if (y == 0)
		t1 = (TTree *)f1->Get("VertexCompositeNtuple;1");
	else
	{
		t1 = (TTree *)dir1->Get("VertexCompositeNtuple;1");
	}
	if (!t1)
	{
		cout << "TTree VertexCompositeNtuple;1 cannot be opened! " << endl;
	}
}
void MC_18::SetupBranches(Int_t x)
{
	if (x == 1)
	{
		t1->SetBranchStatus("*", 0);
		t1->SetBranchStatus("pT", 1);
		t1->SetBranchStatus("candSize", 1);
		t1->SetBranchStatus("centrality", 1);
		t1->SetBranchStatus("eta", 1);
		t1->SetBranchStatus("phi", 1);
		t1->SetBranchStatus("mass", 1);
		t1->SetBranchStatus("pTD1", 1);
		t1->SetBranchStatus("EtaD1", 1);
		t1->SetBranchStatus("PhiD1", 1);
		t1->SetBranchStatus("chargeD1", 1);
		t1->SetBranchStatus("pTD2", 1);
		t1->SetBranchStatus("EtaD2", 1);
		t1->SetBranchStatus("PhiD2", 1);
		t1->SetBranchStatus("chargeD2", 1);
		t1->SetBranchStatus("VtxProb", 1);
		t1->SetBranchStatus("y", 1);
		t1->SetBranchStatus("evtSel", 1);
		t1->SetBranchStatus("bestvtxZ", 1);
		t1->SetBranchStatus("HFsumETPlus", 1);
		t1->SetBranchStatus("HFsumETMinus", 1);
		t1->SetBranchStatus("trigMuon1", 1);
		t1->SetBranchStatus("trigMuon2", 1);
		t1->SetBranchStatus("tightMuon1", 1);
		t1->SetBranchStatus("tightMuon2", 1);
		t1->SetBranchStatus("GlbMuon1", 1);
		t1->SetBranchStatus("GlbMuon2", 1);
		t1->SetBranchStatus("PFMuon1", 1);
		t1->SetBranchStatus("PFMuon2", 1);
		t1->SetBranchStatus("GlbTrkChiD1", 1);
		t1->SetBranchStatus("GlbTrkChiD2", 1);
		t1->SetBranchStatus("nMuonHitD1", 1);
		t1->SetBranchStatus("nMuonHitD2", 1);
		t1->SetBranchStatus("nMatchedStationD1", 1);
		t1->SetBranchStatus("nMatchedStationD2", 1);
		t1->SetBranchStatus("nPixelHitD1", 1);
		t1->SetBranchStatus("nPixelHitD2", 1);
		t1->SetBranchStatus("nTrackerLayerD1", 1);
		t1->SetBranchStatus("nTrackerLayerD2", 1);
		t1->SetBranchStatus("muondXYD1", 1);
		t1->SetBranchStatus("muondXYD2", 1);
		t1->SetBranchStatus("muondZD1", 1);
		t1->SetBranchStatus("muondZD2", 1);
		t1->SetBranchStatus("trigHLT", 1);

		t1->SetBranchAddress("pT", pT);
		t1->SetBranchAddress("centrality", &centrality);
		t1->SetBranchAddress("candSize", &candSize);
		t1->SetBranchAddress("eta", eta);
		t1->SetBranchAddress("phi", phi);
		t1->SetBranchAddress("mass", mass);
		t1->SetBranchAddress("VtxProb", VtxProb);
		t1->SetBranchAddress("pTD1", pTD1);
		t1->SetBranchAddress("EtaD1", EtaD1);
		t1->SetBranchAddress("PhiD1", PhiD1);
		t1->SetBranchAddress("chargeD1", chargeD1);
		t1->SetBranchAddress("pTD2", pTD2);
		t1->SetBranchAddress("EtaD2", EtaD2);
		t1->SetBranchAddress("PhiD2", PhiD2);
		t1->SetBranchAddress("chargeD2", chargeD2);
		t1->SetBranchAddress("y", y);
		t1->SetBranchAddress("evtSel", evtSel);
		t1->SetBranchAddress("bestvtxZ", &bestvtxZ);
		t1->SetBranchAddress("HFsumETPlus", &HFsumETPlus);
		t1->SetBranchAddress("HFsumETMinus", &HFsumETMinus);
		t1->SetBranchAddress("trigMuon1", &trigMuon1);
		t1->SetBranchAddress("trigMuon2", &trigMuon2);
		t1->SetBranchAddress("tightMuon1", tightMuon1);
		t1->SetBranchAddress("tightMuon2", tightMuon2);
		t1->SetBranchAddress("GlbMuon1", GlbMuon1);
		t1->SetBranchAddress("GlbMuon2", GlbMuon2);
		t1->SetBranchAddress("PFMuon1", PFMuon1);
		t1->SetBranchAddress("PFMuon2", PFMuon2);
		t1->SetBranchAddress("GlbTrkChiD1", GlbTrkChiD1);
		t1->SetBranchAddress("GlbTrkChiD2", GlbTrkChiD2);
		t1->SetBranchAddress("nMuonHitD1", nMuonHitD1);
		t1->SetBranchAddress("nMuonHitD2", nMuonHitD2);
		t1->SetBranchAddress("nMatchedStationD1", nMatchedStationD1);
		t1->SetBranchAddress("nMatchedStationD2", nMatchedStationD2);
		t1->SetBranchAddress("nPixelHitD1", nPixelHitD1);
		t1->SetBranchAddress("nPixelHitD2", nPixelHitD2);
		t1->SetBranchAddress("nTrackerLayerD1", nTrackerLayerD1);
		t1->SetBranchAddress("nTrackerLayerD2", nTrackerLayerD2);
		t1->SetBranchAddress("muondXYD1", muondXYD1);
		t1->SetBranchAddress("muondXYD2", muondXYD2);
		t1->SetBranchAddress("muondZD1", muondZD1);
		t1->SetBranchAddress("muondZD2", muondZD2);
		t1->SetBranchAddress("trigHLT", trigHLT);
	}
	else
	{
		// Enable required branches
		t1->SetBranchStatus("*", 0);
		t1->SetBranchStatus("candSize", 1);
		t1->SetBranchStatus("centrality", 1);
		t1->SetBranchStatus("pT", 1);
		t1->SetBranchStatus("eta", 1);
		t1->SetBranchStatus("phi", 1);
		t1->SetBranchStatus("mass", 1);
		t1->SetBranchStatus("VtxProb", 1);
		t1->SetBranchStatus("pTD1", 1);
		t1->SetBranchStatus("EtaD1", 1);
		t1->SetBranchStatus("PhiD1", 1);
		t1->SetBranchStatus("chargeD1", 1);
		t1->SetBranchStatus("pTD2", 1);
		t1->SetBranchStatus("EtaD2", 1);
		t1->SetBranchStatus("PhiD2", 1);
		t1->SetBranchStatus("chargeD2", 1);
		t1->SetBranchStatus("weightLHE_gen", 1);
		t1->SetBranchStatus("candSize_gen", 1);
		t1->SetBranchStatus("weight_gen", 1);
		t1->SetBranchStatus("DecayID_gen", 1);
		t1->SetBranchStatus("pTD1_gen", 1);
		t1->SetBranchStatus("EtaD1_gen", 1);
		t1->SetBranchStatus("PhiD1_gen", 1);
		t1->SetBranchStatus("pTD2_gen", 1);
		t1->SetBranchStatus("EtaD2_gen", 1);
		t1->SetBranchStatus("PhiD2_gen", 1);
		t1->SetBranchStatus("PIDD2_gen", 1);
		t1->SetBranchStatus("PIDD1_gen", 1);
		t1->SetBranchStatus("PIDD2", 1);
		t1->SetBranchStatus("PIDD1", 1);
		t1->SetBranchStatus("PID_gen", 1);
		t1->SetBranchStatus("y", 1);
		t1->SetBranchStatus("evtSel", 1);
		t1->SetBranchStatus("bestvtxZ", 1);
		t1->SetBranchStatus("y_gen", 1);
		t1->SetBranchStatus("pT_gen", 1);
		t1->SetBranchStatus("RecIdx_gen", 1);
		t1->SetBranchStatus("HFsumETPlus", 1);
		t1->SetBranchStatus("HFsumETMinus", 1);
		t1->SetBranchStatus("trigMuon1", 1);
		t1->SetBranchStatus("trigMuon2", 1);
		t1->SetBranchStatus("tightMuon1", 1);
		t1->SetBranchStatus("tightMuon2", 1);
		t1->SetBranchStatus("GlbMuon1", 1);
		t1->SetBranchStatus("GlbMuon2", 1);
		t1->SetBranchStatus("PFMuon1", 1);
		t1->SetBranchStatus("PFMuon2", 1);
		t1->SetBranchStatus("GlbTrkChiD1", 1);
		t1->SetBranchStatus("GlbTrkChiD2", 1);
		t1->SetBranchStatus("nMuonHitD1", 1);
		t1->SetBranchStatus("nMuonHitD2", 1);
		t1->SetBranchStatus("nMatchedStationD1", 1);
		t1->SetBranchStatus("nMatchedStationD2", 1);
		t1->SetBranchStatus("nPixelHitD1", 1);
		t1->SetBranchStatus("nPixelHitD2", 1);
		t1->SetBranchStatus("nTrackerLayerD1", 1);
		t1->SetBranchStatus("nTrackerLayerD2", 1);
		t1->SetBranchStatus("muondXYD1", 1);
		t1->SetBranchStatus("muondXYD2", 1);
		t1->SetBranchStatus("muondZD1", 1);
		t1->SetBranchStatus("muondZD2", 1);
		t1->SetBranchStatus("trigHLT", 1);
		t1->SetBranchStatus("chargeD1_gen", 1);
		t1->SetBranchStatus("chargeD2_gen", 1);
		t1->SetBranchStatus("MotherID_gen",1);
		t1->SetBranchStatus("status_gen",1);
		t1->SetBranchStatus("matchGEN",1);

		// Link to variables
		t1->SetBranchAddress("candSize", &candSize);
		t1->SetBranchAddress("centrality", &centrality);
		t1->SetBranchAddress("pT", pT);
		t1->SetBranchAddress("eta", eta);
		t1->SetBranchAddress("phi", phi);
		t1->SetBranchAddress("mass", mass);
		t1->SetBranchAddress("VtxProb", VtxProb);
		t1->SetBranchAddress("pTD1", pTD1);
		t1->SetBranchAddress("EtaD1", EtaD1);
		t1->SetBranchAddress("PhiD1", PhiD1);
		t1->SetBranchAddress("chargeD1", chargeD1);
		t1->SetBranchAddress("pTD2", pTD2);
		t1->SetBranchAddress("EtaD2", EtaD2);
		t1->SetBranchAddress("PhiD2", PhiD2);
		t1->SetBranchAddress("chargeD2", chargeD2);
		t1->SetBranchAddress("weightLHE_gen", &weightLHE_gen);
		t1->SetBranchAddress("candSize_gen", &candSize_gen);
		t1->SetBranchAddress("weight_gen", &weight_gen);
		t1->SetBranchAddress("DecayID_gen", DecayID_gen);
		t1->SetBranchAddress("pTD1_gen", pTD1_gen);
		t1->SetBranchAddress("EtaD1_gen", EtaD1_gen);
		t1->SetBranchAddress("PhiD1_gen", PhiD1_gen);
		t1->SetBranchAddress("pTD2_gen", pTD2_gen);
		t1->SetBranchAddress("EtaD2_gen", EtaD2_gen);
		t1->SetBranchAddress("PhiD2_gen", PhiD2_gen);
		t1->SetBranchAddress("PIDD2_gen", PIDD2_gen);
		t1->SetBranchAddress("PIDD1_gen", PIDD1_gen);
		t1->SetBranchAddress("PIDD2", PIDD2);
		t1->SetBranchAddress("PIDD1", PIDD1);
		t1->SetBranchAddress("PID_gen", PID_gen);
		t1->SetBranchAddress("y", y);
		t1->SetBranchAddress("evtSel", evtSel);
		t1->SetBranchAddress("bestvtxZ", &bestvtxZ);
		t1->SetBranchAddress("y_gen", y_gen);
		t1->SetBranchAddress("pT_gen", pT_gen);
		t1->SetBranchAddress("RecIdx_gen", RecIdx_gen);
		t1->SetBranchAddress("HFsumETPlus", &HFsumETPlus);
		t1->SetBranchAddress("HFsumETMinus", &HFsumETMinus);
		t1->SetBranchAddress("trigMuon1", &trigMuon1);
		t1->SetBranchAddress("trigMuon2", &trigMuon2);
		t1->SetBranchAddress("tightMuon1", tightMuon1);
		t1->SetBranchAddress("tightMuon2", tightMuon2);
		t1->SetBranchAddress("GlbMuon1", GlbMuon1);
		t1->SetBranchAddress("GlbMuon2", GlbMuon2);
		t1->SetBranchAddress("PFMuon1", PFMuon1);
		t1->SetBranchAddress("PFMuon2", PFMuon2);
		t1->SetBranchAddress("GlbTrkChiD1", GlbTrkChiD1);
		t1->SetBranchAddress("GlbTrkChiD2", GlbTrkChiD2);
		t1->SetBranchAddress("nMuonHitD1", nMuonHitD1);
		t1->SetBranchAddress("nMuonHitD2", nMuonHitD2);
		t1->SetBranchAddress("nMatchedStationD1", nMatchedStationD1);
		t1->SetBranchAddress("nMatchedStationD2", nMatchedStationD2);
		t1->SetBranchAddress("nPixelHitD1", nPixelHitD1);
		t1->SetBranchAddress("nPixelHitD2", nPixelHitD2);
		t1->SetBranchAddress("nTrackerLayerD1", nTrackerLayerD1);
		t1->SetBranchAddress("nTrackerLayerD2", nTrackerLayerD2);
		t1->SetBranchAddress("muondXYD1", muondXYD1);
		t1->SetBranchAddress("muondXYD2", muondXYD2);
		t1->SetBranchAddress("muondZD1", muondZD1);
		t1->SetBranchAddress("muondZD2", muondZD2);
		t1->SetBranchAddress("trigHLT", trigHLT);
		t1->SetBranchAddress("chargeD1_gen", chargeD1_gen);
		t1->SetBranchAddress("chargeD2_gen", chargeD2_gen);
		t1->SetBranchAddress("MotherID_gen",MotherID_gen);
		t1->SetBranchAddress("status_gen",status_gen);
		t1->SetBranchAddress("matchGEN",matchGEN);
	}
}
/*bool MC_18::SetupCut(){
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
}*/
void MC_18::Plot_and_Pull(TH1D *hist, TH1D *hist2, TCanvas *c1, TPad *p1, TPad *p2)
{
	gStyle->SetOptFit(1);
	hist2->GetXaxis()->SetTitle("Invariant Mass(GeV)");
	c1->cd();
	TF1 *fit = new TF1("voigt", "[0]*TMath::Voigt(x-[1]-91.1876,[2],[3]+2.4955)", 60, 120);
	int numBins = hist->GetNbinsX();
	int maxBin = hist->GetMaximumBin();
	fit->SetParameter(0, 60000);
	// fit->SetParameter(1,h_m_1->GetMean()-91.1876);
	fit->SetParameter(1, 0);
	fit->SetParameter(2, hist->GetStdDev());
	fit->SetLineWidth(1);
	fit->SetParameter(3, 2.4955);
	p1->SetMargin(0.1, 0.1, 0, 0.1);
	p2->SetMargin(0.1, 0.1, 0.2, 0);
	p1->Draw();
	p2->Draw();
	p1->cd();
	// p1->SetLogy();
	hist->Fit(fit, "LR");
	hist->Draw();
	Double_t integral = hist->Integral(1, numBins);
	TLatex *text = new TLatex();
	text->SetNDC();
	text->DrawLatex(0.2, 0.8, Form("# of entries is %i", int(integral)));
	p2->cd();
	for (int i = 1; i <= numBins; i++)
	{
		double bincenter = hist->GetXaxis()->GetBinCenter(i);
		double fitresult = fit->Eval(bincenter);
		double binerror = hist->GetBinError(i);
		double bincontent = hist->GetBinContent(i);
		double pullvalue = (bincontent - fitresult) / binerror;
		hist2->SetBinContent(i, pullvalue);
		hist2->SetBinError(i, 1);
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
void MC_18::Plot_and_Pull_gen(TH1D *hist, TH1D *hist2, TCanvas *c1, TPad *p1, TPad *p2)
{
	gStyle->SetOptFit(1);
	hist2->GetXaxis()->SetTitle("Invariant Mass(GeV)");
	c1->cd();
	TF1 *fit = new TF1("voigt", "[0]*TMath::Voigt(x-[1]-91.1876,[2],[3]+2.4955)", 61.1876, 121.1876);
	// TF1 *fit = new TF1("BreitWigner","[0]*TMath::BreitWigner(x, [1]+91.1876, [2]+2.4955)",60,120);
	int numBins = hist->GetNbinsX();
	int maxBin = hist->GetMaximumBin();
	fit->SetParameter(0, 20000);
	// fit->SetParameter(1,h_m_1->GetMean()-91.1876);
	fit->SetParameter(1, 0);
	fit->SetLineWidth(1);
	fit->SetParameter(2, 0);
	p1->SetMargin(0.1, 0.1, 0, 0.1);
	p2->SetMargin(0.1, 0.1, 0.2, 0);
	p1->Draw();
	p2->Draw();
	p1->cd();
	// p1->SetLogy();
	hist->Fit(fit, "LQR");
	hist->Draw();
	Double_t integral = hist->Integral(1, numBins);
	TLatex *text = new TLatex();
	text->SetNDC();
	text->DrawLatex(0.2, 0.8, Form("# of entries is %i", int(integral)));
	p2->cd();
	for (int i = 1; i <= numBins; i++)
	{
		double bincenter = hist->GetXaxis()->GetBinCenter(i);
		double fitresult = fit->Eval(bincenter);
		double binerror = hist->GetBinError(i);
		double bincontent = hist->GetBinContent(i);
		double pullvalue = (bincontent - fitresult) / binerror;
		hist2->SetBinContent(i, pullvalue);
		hist2->SetBinError(i, 1);
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
Float_t MC_18::Calc_Z_gen(int indx)
{
	TLorentzVector muon1;
	Float_t pt_muon_1 = this->pTD1_gen[indx];
	Float_t pt_muon_2 = this->pTD2_gen[indx];
	Float_t eta_muon_1 = this->EtaD1_gen[indx];
	Float_t eta_muon_2 = this->EtaD2_gen[indx];
	Float_t phi_muon_1 = this->PhiD1_gen[indx];
	Float_t phi_muon_2 = this->PhiD2_gen[indx];
	muon1.SetPtEtaPhiM(pt_muon_1, eta_muon_1, phi_muon_1, 0.105);
	TLorentzVector muon2;
	muon2.SetPtEtaPhiM(pt_muon_2, eta_muon_2, phi_muon_2, 0.105);
	TLorentzVector ZBoson = muon1 + muon2;
	Float_t M = ZBoson.M();
	return M;
}

void MC_18::getrapidity(int x)
{
	TLorentzVector v1;
	v1.SetPtEtaPhiM(pT[x], eta[x], phi[x], mass[x]);
	rapidity = v1.Rapidity();
}

int MC_18::getcentrality(int opt)
{
	Double_t hiHF = HFsumETPlus + HFsumETMinus;
	if (opt <= 2 && hiHF >= 5199.95)
		return 0;
	if (opt == 3 && hiHF >= 6242.82)
		return 0;

	Int_t binPos = -1;

	if (opt == 0)
	{
		for (int i = 0; i < 200; ++i)
		{
			if (hiHF >= binTableData[i] && hiHF < binTableData[i + 1])
			{
				binPos = i;
				break;
			}
		}
	}

	if (opt == 1)
	{
		for (int i = 0; i < 200; ++i)
		{
			if (hiHF >= binTableDataUp[i] && hiHF < binTableDataUp[i + 1])
			{
				binPos = i;
				break;
			}
		}
	}

	if (opt == 2)
	{
		for (int i = 0; i < 200; ++i)
		{
			if (hiHF >= binTableDataDown[i] && hiHF < binTableDataDown[i + 1])
			{
				binPos = i;
				break;
			}
		}
	}

	if (opt == 3)
	{
		for (int i = 0; i < 200; i++)
		{
			if (hiHF >= binTableMC_Drum5F[i] && hiHF < binTableMC_Drum5F[i + 1])
			{
				binPos = i;
				break;
			}
		}
	}

	binPos = 200 - 1 - binPos;
	return (Int_t)(200 * ((Double_t)binPos) / (200));
}

void MC_18::forceConsistency(TH2D *top, TH2D *bot)
{
	int nbinsxtop = top->GetNbinsX();
	int nbinsytop = top->GetNbinsY();

	/*for (int i=1; i<=nbinsxtop; i++){
		for (int j=1; j<=nbinsytop; j++){
			if (top->GetBinContent(i,j) > bot->GetBinContent(i,j)){
				cout << "Top > Bot for bin " << i << " " << j << "Set the top == bot" << endl;
				top->SetBinContent(i,j,bot->GetBinContent(i,j));
			}
		}
	}*/
	for (int i = 0; i < top->GetSize(); i++)
	{
		if (top->GetBinContent(i) > bot->GetBinContent(i))
		{
			std::cout << "Warning, fixing an inconsistent bin in: " << top->GetTitle() << ", BinCenter: " << top->GetXaxis()->GetBinCenter(i) << ": " << top->GetBinContent(i) / bot->GetBinContent(i) << std::endl;
			top->SetBinContent(i, bot->GetBinContent(i));
		}
	}
}
Bool_t MC_18::tightMuon1Cut(int x, string type)
{
	if (type == "POG")
	{
		return (GlbMuon1[x] && PFMuon1[x] && (GlbTrkChiD1[x] < 10.) &&
				(nMuonHitD1[x] > 0) && (nMatchedStationD1[x] > 1) &&
				(nPixelHitD1[x] > 0) && (nTrackerLayerD1[x] > 5) &&
				(fabs(muondXYD1[x]) < 0.2) && (fabs(muondZD1[x]) < 0.5));
	}
	else
	{
		std::cout << "[ERROR] Tight MuonID is not defined for " << type << std::endl;
		return false;
	}
}

Bool_t MC_18::tightMuon2Cut(int x, string type)
{
	if (type == "POG")
	{
		return (GlbMuon2[x] && PFMuon2[x] && (GlbTrkChiD2[x] < 10.) &&
				(nMuonHitD2[x] > 0) && (nMatchedStationD2[x] > 1) &&
				(nPixelHitD2[x] > 0) && (nTrackerLayerD2[x] > 5) &&
				(fabs(muondXYD2[x]) < 0.2) && (fabs(muondZD2[x]) < 0.5));
	}
	else
	{
		std::cout << "[ERROR] Tight MuonID is not defined for " << type << std::endl;
		return false;
	}
}
double MC_18::getEfficiency(TEfficiency *e, double y, double pt)
{
	double originalPt = pt;
	if (pt >= 200)
	{
		// std::cout << "Very High-Pt Z (>200) detected. Pt = " << originalPt << " I am pretending it's pt is 199.9 for efficiency purposes!" << std::endl;
		pt = 199.9;
	}

	int bin = e->FindFixBin(y, pt);
	float efficiency = e->GetEfficiency(bin);

	if (efficiency > 0 && efficiency <= 1)
		return efficiency;
	else
	{
		std::cout << "efficiency not in the range [0,1], returning 1!" << std::endl;
		std::cout << "Rapidity: " << y << " Pt: " << originalPt << std::endl;
		return 1;
	}
}
RooFormulaVar *MC_18::CreateRatio(const RooAbsPdf *pdf1, const RooAbsPdf *pdf2, RooRealVar *x)
{
	// Ensure pdf1 and pdf2 share the same variable x
	RooArgList pdfs;
	pdfs.add(*pdf1);
	pdfs.add(*pdf2);

	// Create the ratio formula
	RooFormulaVar *ratio = new RooFormulaVar("ratio", "Ratio PDF", "@0/@1", pdfs);
	return ratio;
}
