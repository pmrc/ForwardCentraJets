// Pedro Cipriano, Mar 2012
// DESY, CMS
// Last Update: 20 Nov 2012
//
// merge_data(output_dir,detail,disp_errors)
// merge the different data datasets in one root file and plots the control distributions

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

#include "common_methods.h"

void show_error_variation2(double *final_errors)
{
//print the error variation on the screen

    cout<<" "<<endl;
    cout<<setiosflags(ios::fixed)<<setprecision(3)<<setfill('0'); //setprecision(3);
    cout<<"Merging unfolded data errors"<<endl;
    cout<<"Observable                  Merged  ->  PYTHIA 6 JetMETTau_2010A | PYTHIA 6 - JetMET_2010A | PYTHIA 8 - JetMETTau_2010A | PYTHIA 8 - JetMET_2010A"<<endl;
    cout<<"Delta Phi                   "<<setw(6)<<final_errors[0]<<"  ->  "<<setw(6)<<final_errors[1]<<"                   | "<<setw(6)<<final_errors[2]<<"                  | "<<setw(6)<<final_errors[3]<<"                     | "<<setw(6)<<final_errors[4]<<endl;
    cout<<"Delta Phi Gap               "<<setw(6)<<final_errors[5]<<"  ->  "<<setw(6)<<final_errors[6]<<"                   | "<<setw(6)<<final_errors[7]<<"                  | "<<setw(6)<<final_errors[8]<<"                     | "<<setw(6)<<final_errors[9]<<endl;
    cout<<"Delta Phi Nogap             "<<setw(6)<<final_errors[10]<<"  ->  "<<setw(6)<<final_errors[11]<<"                   | "<<setw(6)<<final_errors[12]<<"                  | "<<setw(6)<<final_errors[13]<<"                     | "<<setw(6)<<final_errors[14]<<endl;
    cout<<"Leading pT inside gap       "<<setw(6)<<final_errors[15]<<"  ->  "<<setw(6)<<final_errors[16]<<"                   | "<<setw(6)<<final_errors[17]<<"                  | "<<setw(6)<<final_errors[18]<<"                     | "<<setw(6)<<final_errors[19]<<endl;
    cout<<"Leading pT outside gap      "<<setw(6)<<final_errors[20]<<"  ->  "<<setw(6)<<final_errors[21]<<"                   | "<<setw(6)<<final_errors[22]<<"                  | "<<setw(6)<<final_errors[23]<<"                     | "<<setw(6)<<final_errors[24]<<endl;
    cout<<"Delta Phi Deta1             "<<setw(6)<<final_errors[25]<<"  ->  "<<setw(6)<<final_errors[26]<<"                   | "<<setw(6)<<final_errors[27]<<"                  | "<<setw(6)<<final_errors[28]<<"                     | "<<setw(6)<<final_errors[29]<<endl;
    cout<<"Delta Phi Deta2             "<<setw(6)<<final_errors[30]<<"  ->  "<<setw(6)<<final_errors[31]<<"                   | "<<setw(6)<<final_errors[32]<<"                  | "<<setw(6)<<final_errors[33]<<"                     | "<<setw(6)<<final_errors[34]<<endl;
    cout<<"Delta Phi Deta3             "<<setw(6)<<final_errors[35]<<"  ->  "<<setw(6)<<final_errors[36]<<"                   | "<<setw(6)<<final_errors[37]<<"                  | "<<setw(6)<<final_errors[38]<<"                     | "<<setw(6)<<final_errors[39]<<endl;
    cout<<"Delta Phi Deta4             "<<setw(6)<<final_errors[40]<<"  ->  "<<setw(6)<<final_errors[41]<<"                   | "<<setw(6)<<final_errors[42]<<"                  | "<<setw(6)<<final_errors[43]<<"                     | "<<setw(6)<<final_errors[44]<<endl;
    cout<<"Delta Phi Gap Deta1         "<<setw(6)<<final_errors[45]<<"  ->  "<<setw(6)<<final_errors[46]<<"                   | "<<setw(6)<<final_errors[47]<<"                  | "<<setw(6)<<final_errors[48]<<"                     | "<<setw(6)<<final_errors[49]<<endl;
    cout<<"Delta Phi Gap Deta2         "<<setw(6)<<final_errors[50]<<"  ->  "<<setw(6)<<final_errors[51]<<"                   | "<<setw(6)<<final_errors[52]<<"                  | "<<setw(6)<<final_errors[53]<<"                     | "<<setw(6)<<final_errors[54]<<endl;
    cout<<"Delta Phi Gap Deta3         "<<setw(6)<<final_errors[55]<<"  ->  "<<setw(6)<<final_errors[56]<<"                   | "<<setw(6)<<final_errors[57]<<"                  | "<<setw(6)<<final_errors[58]<<"                     | "<<setw(6)<<final_errors[59]<<endl;
    cout<<"Delta Phi Gap Deta4         "<<setw(6)<<final_errors[60]<<"  ->  "<<setw(6)<<final_errors[61]<<"                   | "<<setw(6)<<final_errors[62]<<"                  | "<<setw(6)<<final_errors[63]<<"                     | "<<setw(6)<<final_errors[64]<<endl;
    cout<<"Delta Phi noGap Deta1       "<<setw(6)<<final_errors[65]<<"  ->  "<<setw(6)<<final_errors[66]<<"                   | "<<setw(6)<<final_errors[67]<<"                  | "<<setw(6)<<final_errors[68]<<"                     | "<<setw(6)<<final_errors[69]<<endl;
    cout<<"Delta Phi noGap Deta2       "<<setw(6)<<final_errors[70]<<"  ->  "<<setw(6)<<final_errors[71]<<"                   | "<<setw(6)<<final_errors[72]<<"                  | "<<setw(6)<<final_errors[73]<<"                     | "<<setw(6)<<final_errors[74]<<endl;
    cout<<"Delta Phi noGap Deta3       "<<setw(6)<<final_errors[75]<<"  ->  "<<setw(6)<<final_errors[76]<<"                   | "<<setw(6)<<final_errors[77]<<"                  | "<<setw(6)<<final_errors[78]<<"                     | "<<setw(6)<<final_errors[79]<<endl;
    cout<<"Delta Phi noGap Deta4       "<<setw(6)<<final_errors[80]<<"  ->  "<<setw(6)<<final_errors[81]<<"                   | "<<setw(6)<<final_errors[82]<<"                  | "<<setw(6)<<final_errors[83]<<"                     | "<<setw(6)<<final_errors[84]<<endl;
    cout<<"Leading eta* inside gap     "<<setw(6)<<final_errors[85]<<"  ->  "<<setw(6)<<final_errors[86]<<"                   | "<<setw(6)<<final_errors[87]<<"                  | "<<setw(6)<<final_errors[88]<<"                     | "<<setw(6)<<final_errors[89]<<endl;
    cout<<"Delta eta outside gap       "<<setw(6)<<final_errors[90]<<"  ->  "<<setw(6)<<final_errors[91]<<"                   | "<<setw(6)<<final_errors[92]<<"                  | "<<setw(6)<<final_errors[93]<<"                     | "<<setw(6)<<final_errors[94]<<endl;
    cout<<"Delta Phi Norm              "<<setw(6)<<final_errors[95]<<"  ->  "<<setw(6)<<final_errors[96]<<"                   | "<<setw(6)<<final_errors[97]<<"                  | "<<setw(6)<<final_errors[98]<<"                     | "<<setw(6)<<final_errors[99]<<endl;
    cout<<"Delta Phi Gap Norm          "<<setw(6)<<final_errors[100]<<"  ->  "<<setw(6)<<final_errors[101]<<"                   | "<<setw(6)<<final_errors[102]<<"                  | "<<setw(6)<<final_errors[103]<<"                     | "<<setw(6)<<final_errors[104]<<endl;
    cout<<"Delta Phi NoGap Norm        "<<setw(6)<<final_errors[105]<<"  ->  "<<setw(6)<<final_errors[106]<<"                   | "<<setw(6)<<final_errors[107]<<"                  | "<<setw(6)<<final_errors[108]<<"                     | "<<setw(6)<<final_errors[109]<<endl;
    cout<<"Leading pT Inside Gap Norm  "<<setw(6)<<final_errors[110]<<"  ->  "<<setw(6)<<final_errors[111]<<"                   | "<<setw(6)<<final_errors[112]<<"                  | "<<setw(6)<<final_errors[113]<<"                     | "<<setw(6)<<final_errors[114]<<endl;
    cout<<"Leading pT Outside Gap Norm "<<setw(6)<<final_errors[115]<<"  ->  "<<setw(6)<<final_errors[116]<<"                   | "<<setw(6)<<final_errors[117]<<"                  | "<<setw(6)<<final_errors[118]<<"                     | "<<setw(6)<<final_errors[119]<<endl;
    cout<<"Delta Phi Deta1 Norm        "<<setw(6)<<final_errors[120]<<"  ->  "<<setw(6)<<final_errors[121]<<"                   | "<<setw(6)<<final_errors[122]<<"                  | "<<setw(6)<<final_errors[123]<<"                     | "<<setw(6)<<final_errors[124]<<endl;
    cout<<"Delta Phi Deta2 Norm        "<<setw(6)<<final_errors[125]<<"  ->  "<<setw(6)<<final_errors[126]<<"                   | "<<setw(6)<<final_errors[127]<<"                  | "<<setw(6)<<final_errors[128]<<"                     | "<<setw(6)<<final_errors[129]<<endl;
    cout<<"Delta Phi Deta3 Norm        "<<setw(6)<<final_errors[130]<<"  ->  "<<setw(6)<<final_errors[131]<<"                   | "<<setw(6)<<final_errors[132]<<"                  | "<<setw(6)<<final_errors[133]<<"                     | "<<setw(6)<<final_errors[134]<<endl;
    cout<<"Delta Phi Deta4 Norm        "<<setw(6)<<final_errors[135]<<"  ->  "<<setw(6)<<final_errors[136]<<"                   | "<<setw(6)<<final_errors[137]<<"                  | "<<setw(6)<<final_errors[138]<<"                     | "<<setw(6)<<final_errors[139]<<endl;
    cout<<"Delta Phi Deta1 Gap Norm    "<<setw(6)<<final_errors[140]<<"  ->  "<<setw(6)<<final_errors[141]<<"                   | "<<setw(6)<<final_errors[142]<<"                  | "<<setw(6)<<final_errors[143]<<"                     | "<<setw(6)<<final_errors[144]<<endl;
    cout<<"Delta Phi Deta2 Gap Norm    "<<setw(6)<<final_errors[145]<<"  ->  "<<setw(6)<<final_errors[146]<<"                   | "<<setw(6)<<final_errors[147]<<"                  | "<<setw(6)<<final_errors[148]<<"                     | "<<setw(6)<<final_errors[149]<<endl;
    cout<<"Delta Phi Deta3 Gap Norm    "<<setw(6)<<final_errors[150]<<"  ->  "<<setw(6)<<final_errors[151]<<"                   | "<<setw(6)<<final_errors[152]<<"                  | "<<setw(6)<<final_errors[153]<<"                     | "<<setw(6)<<final_errors[154]<<endl;
    cout<<"Delta Phi Deta4 Gap Norm    "<<setw(6)<<final_errors[155]<<"  ->  "<<setw(6)<<final_errors[156]<<"                   | "<<setw(6)<<final_errors[157]<<"                  | "<<setw(6)<<final_errors[158]<<"                     | "<<setw(6)<<final_errors[159]<<endl;
    cout<<"Delta Phi Deta1 NoGap Norm  "<<setw(6)<<final_errors[160]<<"  ->  "<<setw(6)<<final_errors[161]<<"                   | "<<setw(6)<<final_errors[162]<<"                  | "<<setw(6)<<final_errors[163]<<"                     | "<<setw(6)<<final_errors[164]<<endl;
    cout<<"Delta Phi Deta2 NoGap Norm  "<<setw(6)<<final_errors[165]<<"  ->  "<<setw(6)<<final_errors[166]<<"                   | "<<setw(6)<<final_errors[167]<<"                  | "<<setw(6)<<final_errors[168]<<"                     | "<<setw(6)<<final_errors[169]<<endl;
    cout<<"Delta Phi Deta3 NoGap Norm  "<<setw(6)<<final_errors[170]<<"  ->  "<<setw(6)<<final_errors[171]<<"                   | "<<setw(6)<<final_errors[172]<<"                  | "<<setw(6)<<final_errors[173]<<"                     | "<<setw(6)<<final_errors[174]<<endl;
    cout<<"Delta Phi Deta4 NoGap Norm  "<<setw(6)<<final_errors[175]<<"  ->  "<<setw(6)<<final_errors[176]<<"                   | "<<setw(6)<<final_errors[177]<<"                  | "<<setw(6)<<final_errors[178]<<"                     | "<<setw(6)<<final_errors[179]<<endl;
    cout<<"Leading Eta* Inside Norm    "<<setw(6)<<final_errors[180]<<"  ->  "<<setw(6)<<final_errors[181]<<"                   | "<<setw(6)<<final_errors[182]<<"                  | "<<setw(6)<<final_errors[183]<<"                     | "<<setw(6)<<final_errors[184]<<endl;
    cout<<"Delta Eta Outside Norm      "<<setw(6)<<final_errors[185]<<"  ->  "<<setw(6)<<final_errors[186]<<"                   | "<<setw(6)<<final_errors[187]<<"                  | "<<setw(6)<<final_errors[188]<<"                     | "<<setw(6)<<final_errors[189]<<endl;
}

void plot_histogram(TH1D *merged, TH1D *p6_jetmettau, TH1D *p6_jetmet, TH1D *p8_jetmettau, TH1D *p8_jetmet, string path, string fileout, string legend_position = "top_left", bool logscale = true, bool detail = false)
{
//plots the different datasets and the merged result

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
//declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    if (logscale) { gPad->SetLogy(); }
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    double min = 0.0;
    double max = merged->GetMaximum();
    if (merged->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(merged,detail);
    }
    else
    {
    min = merged->GetMinimum();
    }
    
    set_histogram_min_max(p6_jetmettau, min, max, detail);
    set_histogram_min_max(p6_jetmet, min, max, detail);
    set_histogram_min_max(p8_jetmettau, min, max, detail);
    set_histogram_min_max(p8_jetmet, min, max, detail);
    
    if (logscale)
	{
	max = 1.3 * max;
    	min = 0.7 * min;
    	}
    else
	{
	max = 0.8 * max;
    	min = 0;
	}
    
//plooting
    merged->SetMaximum(max);
    merged->SetMinimum(min);
    merged->SetFillColor(5);
    merged->Draw("e2");
    p6_jetmettau->SetLineWidth(4);
    p6_jetmettau->SetLineColor(2);
    p6_jetmettau->SetLineStyle(1);
    p6_jetmettau->Draw("e1same");
    p6_jetmet->SetLineWidth(4);
    p6_jetmet->SetLineColor(1);
    p6_jetmet->SetLineStyle(2);
    p6_jetmet->Draw("e1same");
    p8_jetmettau->SetLineWidth(4);
    p8_jetmettau->SetLineColor(4);
    p8_jetmettau->SetLineStyle(4);
    p8_jetmettau->Draw("e1same");
    p8_jetmet->SetLineWidth(4);
    p8_jetmet->SetLineColor(6);
    p8_jetmet->SetLineStyle(5);
    p8_jetmet->Draw("e1same");
    
//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 5, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(p6_jetmettau,"PYTHIA 6 - JetMETTau_2010A - L = 284 nb^{-1}","l");
    leg01->AddEntry(p6_jetmet,"PYTHIA 6 - JetMET_2010A - L = 2.90 pb^{-1}","l");
    leg01->AddEntry(p8_jetmettau,"PYTHIA 8 - JetMETTau_2010A - L = 284 nb^{-1}","l");
    leg01->AddEntry(p8_jetmet,"PYTHIA 8 - JetMET_2010A - L = 2.90 pb^{-1}","l");
    leg01->AddEntry(merged,"Merged Data - L = 3.2 pb^{-1}","f");
  //  leg01->AddEntry(merged,"Merged Data - L = 35.2 pb^{-1}","f");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}

void calc_nomalization(TH1D *jetmettau, TH1D *jetmet, double *weigth, bool detail)
{
//computes the merging factors

//declares the variables needed for the calculation
    double entries[2] = {0.,0.};
    double total_entries = 0.0;
    double eff_entries[2] = {0.,0.};
    double eff_total_entries = 0.;
    double integral[2] = {0.,0.};
    
//reset the weigth and scale arrays
    weigth[0] = 1.0;
    weigth[1] = 1.0;

//calculate the scale for the first dataset
    entries[0] = jetmettau->GetEntries();
    if (entries[0] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[0] << endl; } entries[0] = 1.0; }
    eff_entries[0] = jetmettau->GetEffectiveEntries();
    integral[0] = jetmettau->Integral();
    if (integral[0] <= 0) { if (detail) { cout<<"Bad integral = " << integral[0] << endl; } integral[0] = 1.0; }
    if (detail) { cout<<"JetMETTau_2010A "<<setw(8)<<entries[0]<<endl; }

//calculate the scale for the second dataset
    entries[1] = jetmet->GetEntries();
    if (entries[1] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[1] << endl; } entries[1] = 1.0; }
    eff_entries[1] = jetmet->GetEffectiveEntries();
    integral[1] = jetmet->Integral();
    if (integral[1] <= 0) { if (detail) { cout<<"Bad integral = " << integral[1] << endl; } integral[1] = 1.0; }
    if (detail) { cout<<"JetMET_2010A    "<<setw(8)<<entries[1]<<endl; }

//calculate the total number of entries
    total_entries = entries[0] + entries[1];
    eff_total_entries = eff_entries[0] + eff_entries[1];
    if (total_entries <= 0) { if (detail) { cout<<"Bad number of total entries = " << total_entries << endl; } total_entries = 1.0; }

//a few control outputs
    if (detail)
	{
	cout<<"Integral   "<<setw(6)<<integral[0]<<"    "<<setw(6)<<integral[1]<<endl;
	cout<<"Entries    "<<setw(6)<<entries[0]<<"     "<<setw(6)<<entries[1]<<"  "<<setw(6)<<total_entries<<endl;
	}

//calculate the relative weight    
    weigth[0] = entries[0]/total_entries;
    weigth[1] = entries[1]/total_entries;
    if (detail) { cout<<"Weigths "<<weigth[0]<<" "<<weigth[1]<<endl; }

}

void calc_error(TH1D *merged, TH1D *p6_jetmettau, TH1D *p6_jetmet, TH1D *p8_jetmettau, TH1D *p8_jetmet, double *final_errors, int index, bool detail)
{
//evaluate the statistical errors on the merged samples and the individual samples

if (detail) { cout<<"Evaluating statistical errors ..."<<endl; }
//variable declaration
double error_merged = 0.0;
double error[4] = {0.,0.,0.,0.};

//setting the number of bins
int tot_merged = 	merged->GetNbinsX();
int tot_p6_jetmettau = 	p6_jetmettau->GetNbinsX();
int tot_p6_jetmet = 	p6_jetmet->GetNbinsX();
int tot_p8_jetmettau = 	p8_jetmettau->GetNbinsX();
int tot_p8_jetmet = 	p8_jetmet->GetNbinsX();

//remove the empty bins from the bin count
for (int i=1; i<= merged->GetNbinsX();i++)
    {
	if (merged->GetBinContent(i) == 0) 		{ tot_merged = tot_merged - 1; }
	if (p6_jetmettau->GetBinContent(i) == 0) 	{ tot_p6_jetmettau = tot_p6_jetmettau - 1; }
	if (p6_jetmet->GetBinContent(i) == 0) 		{ tot_p6_jetmet = tot_p6_jetmet - 1; }
	if (p8_jetmettau->GetBinContent(i) == 0) 	{ tot_p8_jetmettau = tot_p8_jetmettau - 1; }
	if (p8_jetmet->GetBinContent(i) == 0) 		{ tot_p8_jetmet = tot_p8_jetmet - 1; }
    }

//loop over the bins to get and average the statistical error
    for (int i=1; i<= merged->GetNbinsX();i++)
    {
	if (merged->GetBinContent(i) > 0) { error_merged = error_merged + (merged->GetBinError(i)/merged->GetBinContent(i))/tot_merged; }
	if (p6_jetmettau->GetBinContent(i) > 0) { error[0] = error[0] + (p6_jetmettau->GetBinError(i)/p6_jetmettau->GetBinContent(i))/tot_p6_jetmettau; }
	if (p6_jetmet->GetBinContent(i) > 0) { error[1] = error[1] + (p6_jetmet->GetBinError(i)/p6_jetmet->GetBinContent(i))/tot_p6_jetmet; }
	if (p8_jetmettau->GetBinContent(i) > 0) { error[2] = error[2] + (p8_jetmettau->GetBinError(i)/p8_jetmettau->GetBinContent(i))/tot_p8_jetmettau; }
	if (p8_jetmet->GetBinContent(i) > 0) { error[3] = error[3] + (p8_jetmet->GetBinError(i)/p8_jetmet->GetBinContent(i))/tot_p8_jetmet; }
    }
    if (detail) { cout<<"Errors          "<<error[0]<<" "<<error[1]<<" "<<error[2]<<" "<<error[3]<<" -> "<<error_merged<<endl; }

//saving the errors on final_errors array
    final_errors[index*5+0] = error_merged * 100;
    final_errors[index*5+1] = error[0] * 100;
    final_errors[index*5+2] = error[1] * 100;
    final_errors[index*5+3] = error[2] * 100;
    final_errors[index*5+4] = error[3] * 100;
}


void merge(TH1D *merged, TH1D *p6_jetmettau, TH1D *p6_jetmet, TH1D *p8_jetmettau, TH1D *p8_jetmet, double weigth[2])
{
/*
//try to do a better estimate of correlared errors
//declaring temporary histogram
    TH1D *temp1;
    TH1D *temp2;

//cloning temporary histograms
    temp1 = (TH1D*) merged->Clone();
    temp2 = (TH1D*) merged->Clone();

//merging the different MCs
    temp1->Add(p6_jetmettau,0.5,"B");
    temp1->Add(p8_jetmettau,0.5,"B");
    temp2->Add(p6_jetmet,0.5,"B");
    temp2->Add(p8_jetmet,0.5,"B");

//merges the different datasets
    merged->Add(temp1,weigth[0]);
    merged->Add(temp2,weigth[1]);
*/

    merged->Add(p6_jetmettau,0.5*weigth[0]);
    merged->Add(p8_jetmettau,0.5*weigth[0]);
    merged->Add(p6_jetmet,0.5*weigth[1]);
    merged->Add(p8_jetmet,0.5*weigth[1]);

}

void merge_unfolded_data(string path_p6_jetmettau, string path_p6_jetmet, string path_p8_jetmettau, string path_p8_jetmet, string output_path_plots = "../output/merge_data/", string output_rootfile = "../output/histograms/xsec_data_2010.root", bool detail = false, bool disp_errors = true, bool test = false)
{
//main merging routine
//output detail: true-output all the major steps; false-run quietly
//display errors: true-display merging final errors; false-dont show anything

   if (detail) { cout << "Merge Data Configuration"<<endl; }
   if (detail) { cout << "Input Path PYTHIA 6 JetMETTau_2010A:  "<<path_p6_jetmettau<<endl; }
   if (detail) { cout << "Input Path PYTHIA 6 JetMET_2010A:     "<<path_p6_jetmet<<endl; }
   if (detail) { cout << "Input Path PYTHIA 8 JetMETTau_2010A:  "<<path_p8_jetmettau<<endl; }
   if (detail) { cout << "Input Path PYTHIA 8 JetMET_2010A:     "<<path_p8_jetmet<<endl; }
   if (detail) { cout << "Plot Output Path :                    "<<output_path_plots<<endl; }
   if (detail) { cout << "Root Output Path :                    "<<output_rootfile<<endl; }
   if (detail) { cout << "Detail :                              "<<detail<<endl; }
   if (detail) { cout << "Disp Errors :                         "<<disp_errors<<endl; }
   if (detail) { cout << "Test Mode :                           "<<test<<endl; }

//opening the input and output data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *data_p6_jetmettau 	= new TFile( path_p6_jetmettau.c_str() );
    TFile *data_p6_jetmet 	= new TFile( path_p6_jetmet.c_str() );
    TFile *data_p8_jetmettau 	= new TFile( path_p8_jetmettau.c_str() );
    TFile *data_p8_jetmet 	= new TFile( path_p8_jetmet.c_str() );

//setting the labels for the ploting
    TString label1 = "JetMETTau2010A - L = 284 nb^{-1}";
    TString label2 = "JetMET2010A - L = 2.90 pb^{-1}";
    TString label3 = "Jet2010B - L = 32.1 pb^{-1}";
    TString label4 = "Merged Data - L = 3.18 pb^{-1}";

//declare histogram bins
//int cent_nbins = 7;
//int forw_nbins = 7;
//int all_nbins = 11;

//double cent_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
//double forw_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
//double all_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int in_nbins = 9;
int out_nbins = 9;

double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

//int deta_nbins = 4;
int dphi_nbins = 7;

//double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

//int etac_nbins = 6;
//double etac_bins[7] = {-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8};

//int etaf_nbins = 7;
//double etaf_bins[8] = {-4.7,-4.2,-3.7,-3.2,3.2,3.7,4.2,4.7};

//int eta_nbins = 14;
//double eta_bins[15] = {-4.7,-4.2,-3.7,-3.2,-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8,3.2,3.7,4.2,4.7};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//declaring and setting variables
double weigth[2] = {0.,0.};
double final_errors[38*5];

    for (int i=0; i<= 38*4-1;i++)
    {
    final_errors[i] = 0.0;
    }

//start the merging of the data

//merging the delta phi distribution
    TH1D *delta_phi_p6_jetmettau  = 0;
    TH1D *delta_phi_p6_jetmet  = 0;
    TH1D *delta_phi_p8_jetmettau  = 0;
    TH1D *delta_phi_p8_jetmet  = 0;
    TString delta_phi_name = "output_true_delta_phi";

    data_p6_jetmettau->GetObject(delta_phi_name,delta_phi_p6_jetmettau);
    if (delta_phi_p6_jetmettau == 0) { cout << delta_phi_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_name,delta_phi_p6_jetmet);
    if (delta_phi_p6_jetmet == 0) { cout << delta_phi_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_name,delta_phi_p8_jetmettau);
    if (delta_phi_p8_jetmettau == 0) { cout << delta_phi_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_name,delta_phi_p8_jetmet);
    if (delta_phi_p8_jetmet == 0) { cout << delta_phi_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    TH1D *delta_phi;
    delta_phi = new TH1D("ak5PF_delta_phi","#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi->Sumw2();

    if (detail) { cout<<"ak5PF_delta Phi"<<endl; }
    calc_nomalization(delta_phi_p6_jetmettau, delta_phi_p6_jetmet, weigth, detail);
    merge(delta_phi, delta_phi_p6_jetmettau, delta_phi_p6_jetmet, delta_phi_p8_jetmettau, delta_phi_p8_jetmet, weigth);    
    plot_histogram(delta_phi, delta_phi_p6_jetmettau, delta_phi_p6_jetmet, delta_phi_p8_jetmettau, delta_phi_p8_jetmet, output_path_plots, "merging_delta_phi", "top_left", true, detail);
    plot_histogram(delta_phi, delta_phi_p6_jetmettau, delta_phi_p6_jetmet, delta_phi_p8_jetmettau, delta_phi_p8_jetmet, output_path_plots, "merging_delta_phi_linear", "top_left", false, detail);
    calc_error(delta_phi, delta_phi_p6_jetmettau, delta_phi_p6_jetmet, delta_phi_p8_jetmettau, delta_phi_p8_jetmet, final_errors, 0, detail);    
    

//merging the delta phi norm distribution
    TH1D *delta_phi_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_norm_p6_jetmet  = 0;
    TH1D *delta_phi_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_norm_p8_jetmet  = 0;
    TString delta_phi_norm_name = "output_true_delta_phi_norm";

    data_p6_jetmettau->GetObject(delta_phi_norm_name,delta_phi_norm_p6_jetmettau);
    if (delta_phi_norm_p6_jetmettau == 0) { cout << delta_phi_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_norm_name,delta_phi_norm_p6_jetmet);
    if (delta_phi_norm_p6_jetmet == 0) { cout << delta_phi_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_norm_name,delta_phi_norm_p8_jetmettau);
    if (delta_phi_norm_p8_jetmettau == 0) { cout << delta_phi_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_norm_name,delta_phi_norm_p8_jetmet);
    if (delta_phi_norm_p8_jetmet == 0) { cout << delta_phi_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    TH1D *delta_phi_norm;
    delta_phi_norm = new TH1D("ak5PF_delta_phi_norm","#Delta#phi;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_norm->Sumw2();

    if (detail) { cout<<"ak5PF_delta Phi Norm"<<endl; }
    calc_nomalization(delta_phi_norm_p6_jetmettau, delta_phi_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_norm, delta_phi_norm_p6_jetmettau, delta_phi_norm_p6_jetmet, delta_phi_norm_p8_jetmettau, delta_phi_norm_p8_jetmet, weigth);    
    plot_histogram(delta_phi_norm, delta_phi_norm_p6_jetmettau, delta_phi_norm_p6_jetmet, delta_phi_norm_p8_jetmettau, delta_phi_norm_p8_jetmet, output_path_plots, "merging_delta_phi_norm", "top_left", true, detail);
    plot_histogram(delta_phi_norm, delta_phi_norm_p6_jetmettau, delta_phi_norm_p6_jetmet, delta_phi_norm_p8_jetmettau, delta_phi_norm_p8_jetmet, output_path_plots, "merging_delta_phi_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_norm, delta_phi_norm_p6_jetmettau, delta_phi_norm_p6_jetmet, delta_phi_norm_p8_jetmettau, delta_phi_norm_p8_jetmet, final_errors, 19, detail);   


//merging the delta phi with gap distribution
    TH1D *delta_phi_gap_p6_jetmettau  = 0;
    TH1D *delta_phi_gap_p6_jetmet  = 0;
    TH1D *delta_phi_gap_p8_jetmettau  = 0;
    TH1D *delta_phi_gap_p8_jetmet  = 0;
    TString delta_phi_gap_name = "output_true_delta_phi_gap";

    data_p6_jetmettau->GetObject(delta_phi_gap_name,delta_phi_gap_p6_jetmettau);
    if (delta_phi_gap_p6_jetmettau == 0) { cout << delta_phi_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_gap_name,delta_phi_gap_p6_jetmet);
    if (delta_phi_gap_p6_jetmet == 0) { cout << delta_phi_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_gap_name,delta_phi_gap_p8_jetmettau);
    if (delta_phi_gap_p8_jetmettau == 0) { cout << delta_phi_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_gap_name,delta_phi_gap_p8_jetmet);
    if (delta_phi_gap_p8_jetmet == 0) { cout << delta_phi_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D("ak5PF_delta_phi_gap","#Delta#phi^{gap};|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_gap->Sumw2();

    if (detail) { cout<<"ak5PF_delta Phi gap"<<endl; }
    calc_nomalization(delta_phi_gap_p6_jetmettau, delta_phi_gap_p6_jetmet, weigth, detail);
    merge(delta_phi_gap, delta_phi_gap_p6_jetmettau, delta_phi_gap_p6_jetmet, delta_phi_gap_p8_jetmettau, delta_phi_gap_p8_jetmet, weigth);    
    plot_histogram(delta_phi_gap, delta_phi_gap_p6_jetmettau, delta_phi_gap_p6_jetmet, delta_phi_gap_p8_jetmettau, delta_phi_gap_p8_jetmet, output_path_plots, "merging_delta_phi_gap", "top_left", true, detail);
    plot_histogram(delta_phi_gap, delta_phi_gap_p6_jetmettau, delta_phi_gap_p6_jetmet, delta_phi_gap_p8_jetmettau, delta_phi_gap_p8_jetmet, output_path_plots, "merging_delta_phi_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_gap, delta_phi_gap_p6_jetmettau, delta_phi_gap_p6_jetmet, delta_phi_gap_p8_jetmettau, delta_phi_gap_p8_jetmet, final_errors, 1, detail); 



//merging the delta phi with gap norm distribution
    TH1D *delta_phi_gap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_gap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_gap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_gap_norm_p8_jetmet  = 0;
    TString delta_phi_gap_norm_name = "output_true_delta_phi_gap_norm";

    data_p6_jetmettau->GetObject(delta_phi_gap_norm_name,delta_phi_gap_norm_p6_jetmettau);
    if (delta_phi_gap_norm_p6_jetmettau == 0) { cout << delta_phi_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_gap_norm_name,delta_phi_gap_norm_p6_jetmet);
    if (delta_phi_gap_norm_p6_jetmet == 0) { cout << delta_phi_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_gap_norm_name,delta_phi_gap_norm_p8_jetmettau);
    if (delta_phi_gap_norm_p8_jetmettau == 0) { cout << delta_phi_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_gap_norm_name,delta_phi_gap_norm_p8_jetmet);
    if (delta_phi_gap_norm_p8_jetmet == 0) { cout << delta_phi_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_gap_norm;
    delta_phi_gap_norm =  new TH1D("ak5PF_delta_phi_gap_norm","#Delta#phi^{gap};|#Delta#phi| [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_gap_norm->Sumw2();

    if (detail) { cout<<"ak5PF_delta Phi Gap Norm"<<endl; }
    calc_nomalization(delta_phi_gap_norm_p6_jetmettau, delta_phi_gap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_gap_norm, delta_phi_gap_norm_p6_jetmettau, delta_phi_gap_norm_p6_jetmet, delta_phi_gap_norm_p8_jetmettau, delta_phi_gap_norm_p8_jetmet, weigth);    
    plot_histogram(delta_phi_gap_norm, delta_phi_gap_norm_p6_jetmettau, delta_phi_gap_norm_p6_jetmet, delta_phi_gap_norm_p8_jetmettau, delta_phi_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_gap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_gap_norm, delta_phi_gap_norm_p6_jetmettau, delta_phi_gap_norm_p6_jetmet, delta_phi_gap_norm_p8_jetmettau, delta_phi_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_gap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_gap_norm, delta_phi_gap_norm_p6_jetmettau, delta_phi_gap_norm_p6_jetmet, delta_phi_gap_norm_p8_jetmettau, delta_phi_gap_norm_p8_jetmet, final_errors, 20, detail); 


//merging the delta phi without gap distribution 
    TH1D *delta_phi_nogap_p6_jetmettau  = 0;
    TH1D *delta_phi_nogap_p6_jetmet  = 0;
    TH1D *delta_phi_nogap_p8_jetmettau  = 0;
    TH1D *delta_phi_nogap_p8_jetmet  = 0;
    TString delta_phi_nogap_name = "output_true_delta_phi_nogap";

    data_p6_jetmettau->GetObject(delta_phi_nogap_name,delta_phi_nogap_p6_jetmettau);
    if (delta_phi_nogap_p6_jetmettau == 0) { cout << delta_phi_nogap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_nogap_name,delta_phi_nogap_p6_jetmet);
    if (delta_phi_nogap_p6_jetmet == 0) { cout << delta_phi_nogap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_nogap_name,delta_phi_nogap_p8_jetmettau);
    if (delta_phi_nogap_p8_jetmettau == 0) { cout << delta_phi_nogap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_nogap_name,delta_phi_nogap_p8_jetmet);
    if (delta_phi_nogap_p8_jetmet == 0) { cout << delta_phi_nogap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D("ak5PF_delta_phi_nogap","#Delta#phi^{no gap};|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Phi nogap"<<endl; }
    calc_nomalization(delta_phi_nogap_p6_jetmettau, delta_phi_nogap_p6_jetmet, weigth, detail);
    merge(delta_phi_nogap, delta_phi_nogap_p6_jetmettau, delta_phi_nogap_p6_jetmet, delta_phi_nogap_p8_jetmettau, delta_phi_nogap_p8_jetmet, weigth);    
    plot_histogram(delta_phi_nogap, delta_phi_nogap_p6_jetmettau, delta_phi_nogap_p6_jetmet, delta_phi_nogap_p8_jetmettau, delta_phi_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_nogap, delta_phi_nogap_p6_jetmettau, delta_phi_nogap_p6_jetmet, delta_phi_nogap_p8_jetmettau, delta_phi_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_nogap, delta_phi_nogap_p6_jetmettau, delta_phi_nogap_p6_jetmet, delta_phi_nogap_p8_jetmettau, delta_phi_nogap_p8_jetmet, final_errors, 2, detail); 


//merging the delta phi without gap norm distribution 
    TH1D *delta_phi_nogap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_nogap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_nogap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_nogap_norm_p8_jetmet  = 0;
    TString delta_phi_nogap_norm_name = "output_true_delta_phi_nogap_norm";

    data_p6_jetmettau->GetObject(delta_phi_nogap_norm_name,delta_phi_nogap_norm_p6_jetmettau);
    if (delta_phi_nogap_norm_p6_jetmettau == 0) { cout << delta_phi_nogap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_nogap_norm_name,delta_phi_nogap_norm_p6_jetmet);
    if (delta_phi_nogap_norm_p6_jetmet == 0) { cout << delta_phi_nogap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_nogap_norm_name,delta_phi_nogap_norm_p8_jetmettau);
    if (delta_phi_nogap_norm_p8_jetmettau == 0) { cout << delta_phi_nogap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_nogap_norm_name,delta_phi_nogap_norm_p8_jetmet);
    if (delta_phi_nogap_norm_p8_jetmet == 0) { cout << delta_phi_nogap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_nogap_norm;
    delta_phi_nogap_norm =  new TH1D("ak5PF_delta_phi_nogap_norm","#Delta#phi^{Nogap} [rad];#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_nogap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Phi Nogap Norm"<<endl; }
    calc_nomalization(delta_phi_nogap_norm_p6_jetmettau, delta_phi_nogap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_nogap_norm, delta_phi_nogap_norm_p6_jetmettau, delta_phi_nogap_norm_p6_jetmet, delta_phi_nogap_norm_p8_jetmettau, delta_phi_nogap_norm_p8_jetmet, weigth);    
    plot_histogram(delta_phi_nogap_norm, delta_phi_nogap_norm_p6_jetmettau, delta_phi_nogap_norm_p6_jetmet, delta_phi_nogap_norm_p8_jetmettau, delta_phi_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_nogap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_nogap_norm, delta_phi_nogap_norm_p6_jetmettau, delta_phi_nogap_norm_p6_jetmet, delta_phi_nogap_norm_p8_jetmettau, delta_phi_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_nogap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_nogap_norm, delta_phi_nogap_norm_p6_jetmettau, delta_phi_nogap_norm_p6_jetmet, delta_phi_nogap_norm_p8_jetmettau, delta_phi_nogap_norm_p8_jetmet, final_errors, 21, detail); 


//merging the leading pt inside the gap distribution 
    TH1D *leading_pt_inside_gap_p6_jetmettau  = 0;
    TH1D *leading_pt_inside_gap_p6_jetmet  = 0;
    TH1D *leading_pt_inside_gap_p8_jetmettau  = 0;
    TH1D *leading_pt_inside_gap_p8_jetmet  = 0;
    TString leading_pt_inside_gap_name = "output_true_leading_pt_inside_gap";

    data_p6_jetmettau->GetObject(leading_pt_inside_gap_name,leading_pt_inside_gap_p6_jetmettau);
    if (leading_pt_inside_gap_p6_jetmettau == 0) { cout << leading_pt_inside_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(leading_pt_inside_gap_name,leading_pt_inside_gap_p6_jetmet);
    if (leading_pt_inside_gap_p6_jetmet == 0) { cout << leading_pt_inside_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(leading_pt_inside_gap_name,leading_pt_inside_gap_p8_jetmettau);
    if (leading_pt_inside_gap_p8_jetmettau == 0) { cout << leading_pt_inside_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(leading_pt_inside_gap_name,leading_pt_inside_gap_p8_jetmet);
    if (leading_pt_inside_gap_p8_jetmet == 0) { cout << leading_pt_inside_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap = new TH1D("ak5PF_leading_pt_inside_gap","Leading_pt_inside_gap;p_{T} [GeV];#frac{d#sigma}{dp_{T}}", in_nbins, in_bins);
    leading_pt_inside_gap->Sumw2();
    
    if (detail) { cout<<"Leading pt inside gap"<<endl; }
    calc_nomalization(leading_pt_inside_gap_p6_jetmettau, leading_pt_inside_gap_p6_jetmet, weigth, detail);
    merge(leading_pt_inside_gap, leading_pt_inside_gap_p6_jetmettau, leading_pt_inside_gap_p6_jetmet, leading_pt_inside_gap_p8_jetmettau, leading_pt_inside_gap_p8_jetmet, weigth);
    plot_histogram(leading_pt_inside_gap, leading_pt_inside_gap_p6_jetmettau, leading_pt_inside_gap_p6_jetmet, leading_pt_inside_gap_p8_jetmettau, leading_pt_inside_gap_p8_jetmet, output_path_plots, "merging_leading_pt_inside_gap", "top_right", true, detail);
    plot_histogram(leading_pt_inside_gap, leading_pt_inside_gap_p6_jetmettau, leading_pt_inside_gap_p6_jetmet, leading_pt_inside_gap_p8_jetmettau, leading_pt_inside_gap_p8_jetmet, output_path_plots, "merging_leading_pt_inside_gap_linear", "top_right", false, detail);
    calc_error(leading_pt_inside_gap, leading_pt_inside_gap_p6_jetmettau, leading_pt_inside_gap_p6_jetmet, leading_pt_inside_gap_p8_jetmettau, leading_pt_inside_gap_p8_jetmet, final_errors, 3, detail);


//merging the leading pt inside the gap Norm distribution 
    TH1D *leading_pt_inside_gap_norm_p6_jetmettau  = 0;
    TH1D *leading_pt_inside_gap_norm_p6_jetmet  = 0;
    TH1D *leading_pt_inside_gap_norm_p8_jetmettau  = 0;
    TH1D *leading_pt_inside_gap_norm_p8_jetmet  = 0;
    TString leading_pt_inside_gap_norm_name = "output_true_leading_pt_inside_gap_norm";

    data_p6_jetmettau->GetObject(leading_pt_inside_gap_norm_name,leading_pt_inside_gap_norm_p6_jetmettau);
    if (leading_pt_inside_gap_norm_p6_jetmettau == 0) { cout << leading_pt_inside_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(leading_pt_inside_gap_norm_name,leading_pt_inside_gap_norm_p6_jetmet);
    if (leading_pt_inside_gap_norm_p6_jetmet == 0) { cout << leading_pt_inside_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(leading_pt_inside_gap_norm_name,leading_pt_inside_gap_norm_p8_jetmettau);
    if (leading_pt_inside_gap_norm_p8_jetmettau == 0) { cout << leading_pt_inside_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(leading_pt_inside_gap_norm_name,leading_pt_inside_gap_norm_p8_jetmet);
    if (leading_pt_inside_gap_norm_p8_jetmet == 0) { cout << leading_pt_inside_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    
    TH1D *leading_pt_inside_gap_norm;
    leading_pt_inside_gap_norm = new TH1D("ak5PF_leading_pt_inside_gap_norm","Leading_pt_inside_gap;p_{T}^{inside} [GeV];#frac{d#sigma}{dp_{T}^{inside}}", in_nbins, in_bins);
    leading_pt_inside_gap_norm->Sumw2();
    
    if (detail) { cout<<"Leading pt inside Gap Norm"<<endl; }
    calc_nomalization(leading_pt_inside_gap_norm_p6_jetmettau, leading_pt_inside_gap_norm_p6_jetmet, weigth, detail);
    merge(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_p6_jetmettau, leading_pt_inside_gap_norm_p6_jetmet, leading_pt_inside_gap_norm_p8_jetmettau, leading_pt_inside_gap_norm_p8_jetmet, weigth);
    plot_histogram(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_p6_jetmettau, leading_pt_inside_gap_norm_p6_jetmet, leading_pt_inside_gap_norm_p8_jetmettau, leading_pt_inside_gap_norm_p8_jetmet, output_path_plots, "merging_leading_pt_inside_gap_norm", "top_right", true, detail);
    plot_histogram(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_p6_jetmettau, leading_pt_inside_gap_norm_p6_jetmet, leading_pt_inside_gap_norm_p8_jetmettau, leading_pt_inside_gap_norm_p8_jetmet, output_path_plots, "merging_leading_pt_inside_gap_norm_linear", "top_right", false, detail);
    calc_error(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_p6_jetmettau, leading_pt_inside_gap_norm_p6_jetmet, leading_pt_inside_gap_norm_p8_jetmettau, leading_pt_inside_gap_norm_p8_jetmet, final_errors, 22, detail);


//merging the leading pt outside the gap distribution    
    TH1D *leading_pt_outside_gap_p6_jetmettau  = 0;
    TH1D *leading_pt_outside_gap_p6_jetmet  = 0;
    TH1D *leading_pt_outside_gap_p8_jetmettau  = 0;
    TH1D *leading_pt_outside_gap_p8_jetmet  = 0;
    TString leading_pt_outside_gap_name = "output_true_leading_pt_outside_gap";

    data_p6_jetmettau->GetObject(leading_pt_outside_gap_name,leading_pt_outside_gap_p6_jetmettau);
    if (leading_pt_outside_gap_p6_jetmettau == 0) { cout << leading_pt_outside_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(leading_pt_outside_gap_name,leading_pt_outside_gap_p6_jetmet);
    if (leading_pt_outside_gap_p6_jetmet == 0) { cout << leading_pt_outside_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(leading_pt_outside_gap_name,leading_pt_outside_gap_p8_jetmettau);
    if (leading_pt_outside_gap_p8_jetmettau == 0) { cout << leading_pt_outside_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(leading_pt_outside_gap_name,leading_pt_outside_gap_p8_jetmet);
    if (leading_pt_outside_gap_p8_jetmet == 0) { cout << leading_pt_outside_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap = new TH1D("ak5PF_leading_pt_outside_gap","Leading_pt_outside_gap;p_{T} [GeV];#frac{d#sigma}{dp_{T}}", out_nbins, out_bins);
    leading_pt_outside_gap->Sumw2();
    
    if (detail) { cout<<"Leading pt outside gap "<<endl; }
    calc_nomalization(leading_pt_outside_gap_p6_jetmettau, leading_pt_outside_gap_p6_jetmet, weigth, detail);
    merge(leading_pt_outside_gap, leading_pt_outside_gap_p6_jetmettau, leading_pt_outside_gap_p6_jetmet, leading_pt_outside_gap_p8_jetmettau, leading_pt_outside_gap_p8_jetmet, weigth);
    plot_histogram(leading_pt_outside_gap, leading_pt_outside_gap_p6_jetmettau, leading_pt_outside_gap_p6_jetmet, leading_pt_outside_gap_p8_jetmettau, leading_pt_outside_gap_p8_jetmet, output_path_plots, "merging_leading_pt_outside_gap", "top_right", true, detail);
    plot_histogram(leading_pt_outside_gap, leading_pt_outside_gap_p6_jetmettau, leading_pt_outside_gap_p6_jetmet, leading_pt_outside_gap_p8_jetmettau, leading_pt_outside_gap_p8_jetmet, output_path_plots, "merging_leading_pt_outside_gap_linear", "top_right", false, detail);
    calc_error(leading_pt_outside_gap, leading_pt_outside_gap_p6_jetmettau, leading_pt_outside_gap_p6_jetmet, leading_pt_outside_gap_p8_jetmettau, leading_pt_outside_gap_p8_jetmet, final_errors, 4, detail);


//merging the leading pt outside the gap Norm distribution    
    TH1D *leading_pt_outside_gap_norm_p6_jetmettau  = 0;
    TH1D *leading_pt_outside_gap_norm_p6_jetmet  = 0;
    TH1D *leading_pt_outside_gap_norm_p8_jetmettau  = 0;
    TH1D *leading_pt_outside_gap_norm_p8_jetmet  = 0;
    TString leading_pt_outside_gap_norm_name = "output_true_leading_pt_outside_gap_norm";

    data_p6_jetmettau->GetObject(leading_pt_outside_gap_norm_name,leading_pt_outside_gap_norm_p6_jetmettau);
    if (leading_pt_outside_gap_norm_p6_jetmettau == 0) { cout << leading_pt_outside_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(leading_pt_outside_gap_norm_name,leading_pt_outside_gap_norm_p6_jetmet);
    if (leading_pt_outside_gap_norm_p6_jetmet == 0) { cout << leading_pt_outside_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(leading_pt_outside_gap_norm_name,leading_pt_outside_gap_norm_p8_jetmettau);
    if (leading_pt_outside_gap_norm_p8_jetmettau == 0) { cout << leading_pt_outside_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(leading_pt_outside_gap_norm_name,leading_pt_outside_gap_norm_p8_jetmet);
    if (leading_pt_outside_gap_norm_p8_jetmet == 0) { cout << leading_pt_outside_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *leading_pt_outside_gap_norm;
    leading_pt_outside_gap_norm = new TH1D("ak5PF_leading_pt_outside_gap_norm","Leading_pt_outside_gap;p_{T}^{outside} [GeV];#frac{d#sigma}{dp_{T}^{outside}}", out_nbins, out_bins);
    leading_pt_outside_gap_norm->Sumw2();
    
    if (detail) { cout<<"Leading pt outside Gap Norm"<<endl; }
    calc_nomalization(leading_pt_outside_gap_norm_p6_jetmettau, leading_pt_outside_gap_norm_p6_jetmet, weigth, detail);
    merge(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_p6_jetmettau, leading_pt_outside_gap_norm_p6_jetmet, leading_pt_outside_gap_norm_p8_jetmettau, leading_pt_outside_gap_norm_p8_jetmet, weigth);
    plot_histogram(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_p6_jetmettau, leading_pt_outside_gap_norm_p6_jetmet, leading_pt_outside_gap_norm_p8_jetmettau, leading_pt_outside_gap_norm_p8_jetmet, output_path_plots, "merging_leading_pt_outside_gap_norm", "top_right", true, detail);
    plot_histogram(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_p6_jetmettau, leading_pt_outside_gap_norm_p6_jetmet, leading_pt_outside_gap_norm_p8_jetmettau, leading_pt_outside_gap_norm_p8_jetmet, output_path_plots, "merging_leading_pt_outside_gap_norm_linear", "top_right", false, detail);
    calc_error(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_p6_jetmettau, leading_pt_outside_gap_norm_p6_jetmet, leading_pt_outside_gap_norm_p8_jetmettau, leading_pt_outside_gap_norm_p8_jetmet, final_errors, 23, detail);


//merging the delta phi deta1 distributions
    TH1D *delta_phi_deta1_p6_jetmettau  = 0;
    TH1D *delta_phi_deta1_p6_jetmet  = 0;
    TH1D *delta_phi_deta1_p8_jetmettau  = 0;
    TH1D *delta_phi_deta1_p8_jetmet  = 0;
    TString delta_phi_deta1_name = "output_true_delta_phi_deta1";

    data_p6_jetmettau->GetObject(delta_phi_deta1_name,delta_phi_deta1_p6_jetmettau);
    if (delta_phi_deta1_p6_jetmettau == 0) { cout << delta_phi_deta1_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta1_name,delta_phi_deta1_p6_jetmet);
    if (delta_phi_deta1_p6_jetmet == 0) { cout << delta_phi_deta1_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta1_name,delta_phi_deta1_p8_jetmettau);
    if (delta_phi_deta1_p8_jetmettau == 0) { cout << delta_phi_deta1_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta1_name,delta_phi_deta1_p8_jetmet);
    if (delta_phi_deta1_p8_jetmet == 0) { cout << delta_phi_deta1_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 = new TH1D("ak5PF_delta_phi_deta1","Delta_phi_deta1;|#Delta#phi| [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta1"<<endl; }
    calc_nomalization(delta_phi_deta1_p6_jetmettau, delta_phi_deta1_p6_jetmet, weigth, detail);
    merge(delta_phi_deta1, delta_phi_deta1_p6_jetmettau, delta_phi_deta1_p6_jetmet, delta_phi_deta1_p8_jetmettau, delta_phi_deta1_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta1, delta_phi_deta1_p6_jetmettau, delta_phi_deta1_p6_jetmet, delta_phi_deta1_p8_jetmettau, delta_phi_deta1_p8_jetmet, output_path_plots, "merging_delta_phi_deta1", "top_left", true, detail);
    plot_histogram(delta_phi_deta1, delta_phi_deta1_p6_jetmettau, delta_phi_deta1_p6_jetmet, delta_phi_deta1_p8_jetmettau, delta_phi_deta1_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1, delta_phi_deta1_p6_jetmettau, delta_phi_deta1_p6_jetmet, delta_phi_deta1_p8_jetmettau, delta_phi_deta1_p8_jetmet, final_errors, 5, detail);


//merging the delta phi deta1 norm distributions
    TH1D *delta_phi_deta1_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta1_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta1_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta1_norm_p8_jetmet  = 0;
    TString delta_phi_deta1_norm_name = "output_true_delta_phi_deta1_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta1_norm_name,delta_phi_deta1_norm_p6_jetmettau);
    if (delta_phi_deta1_norm_p6_jetmettau == 0) { cout << delta_phi_deta1_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta1_norm_name,delta_phi_deta1_norm_p6_jetmet);
    if (delta_phi_deta1_norm_p6_jetmet == 0) { cout << delta_phi_deta1_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta1_norm_name,delta_phi_deta1_norm_p8_jetmettau);
    if (delta_phi_deta1_norm_p8_jetmettau == 0) { cout << delta_phi_deta1_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta1_norm_name,delta_phi_deta1_norm_p8_jetmet);
    if (delta_phi_deta1_norm_p8_jetmet == 0) { cout << delta_phi_deta1_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_deta1_norm;
    delta_phi_deta1_norm = new TH1D("ak5PF_delta_phi_deta1_norm","Delta_phi_deta1;#Delta#phi [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma^{2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta1 Norm"<<endl; }
    calc_nomalization(delta_phi_deta1_norm_p6_jetmettau, delta_phi_deta1_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta1_norm, delta_phi_deta1_norm_p6_jetmettau, delta_phi_deta1_norm_p6_jetmet, delta_phi_deta1_norm_p8_jetmettau, delta_phi_deta1_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta1_norm, delta_phi_deta1_norm_p6_jetmettau, delta_phi_deta1_norm_p6_jetmet, delta_phi_deta1_norm_p8_jetmettau, delta_phi_deta1_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_norm, delta_phi_deta1_norm_p6_jetmettau, delta_phi_deta1_norm_p6_jetmet, delta_phi_deta1_norm_p8_jetmettau, delta_phi_deta1_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1_norm, delta_phi_deta1_norm_p6_jetmettau, delta_phi_deta1_norm_p6_jetmet, delta_phi_deta1_norm_p8_jetmettau, delta_phi_deta1_norm_p8_jetmet, final_errors, 24, detail);

//merging the delta phi deta2 distributions    
    TH1D *delta_phi_deta2_p6_jetmettau  = 0;
    TH1D *delta_phi_deta2_p6_jetmet  = 0;
    TH1D *delta_phi_deta2_p8_jetmettau  = 0;
    TH1D *delta_phi_deta2_p8_jetmet  = 0;
    TString delta_phi_deta2_name = "output_true_delta_phi_deta2";

    data_p6_jetmettau->GetObject(delta_phi_deta2_name,delta_phi_deta2_p6_jetmettau);
    if (delta_phi_deta2_p6_jetmettau == 0) { cout << delta_phi_deta2_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta2_name,delta_phi_deta2_p6_jetmet);
    if (delta_phi_deta2_p6_jetmet == 0) { cout << delta_phi_deta2_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta2_name,delta_phi_deta2_p8_jetmettau);
    if (delta_phi_deta2_p8_jetmettau == 0) { cout << delta_phi_deta2_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta2_name,delta_phi_deta2_p8_jetmet);
    if (delta_phi_deta2_p8_jetmet == 0) { cout << delta_phi_deta2_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 = new TH1D("ak5PF_delta_phi_deta2","Delta_phi_deta2;|#Delta#phi| [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta2"<<endl; }
    calc_nomalization(delta_phi_deta2_p6_jetmettau, delta_phi_deta2_p6_jetmet, weigth, detail);
    merge(delta_phi_deta2, delta_phi_deta2_p6_jetmettau, delta_phi_deta2_p6_jetmet, delta_phi_deta2_p8_jetmettau, delta_phi_deta2_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta2, delta_phi_deta2_p6_jetmettau, delta_phi_deta2_p6_jetmet, delta_phi_deta2_p8_jetmettau, delta_phi_deta2_p8_jetmet, output_path_plots, "merging_delta_phi_deta2", "top_left", true, detail);
    plot_histogram(delta_phi_deta2, delta_phi_deta2_p6_jetmettau, delta_phi_deta2_p6_jetmet, delta_phi_deta2_p8_jetmettau, delta_phi_deta2_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2, delta_phi_deta2_p6_jetmettau, delta_phi_deta2_p6_jetmet, delta_phi_deta2_p8_jetmettau, delta_phi_deta2_p8_jetmet, final_errors, 6, detail);

//merging the delta phi deta2 norm distributions    
    TH1D *delta_phi_deta2_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta2_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta2_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta2_norm_p8_jetmet  = 0;
    TString delta_phi_deta2_norm_name = "output_true_delta_phi_deta2_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta2_norm_name,delta_phi_deta2_norm_p6_jetmettau);
    if (delta_phi_deta2_norm_p6_jetmettau == 0) { cout << delta_phi_deta2_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta2_norm_name,delta_phi_deta2_norm_p6_jetmet);
    if (delta_phi_deta2_norm_p6_jetmet == 0) { cout << delta_phi_deta2_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta2_norm_name,delta_phi_deta2_norm_p8_jetmettau);
    if (delta_phi_deta2_norm_p8_jetmettau == 0) { cout << delta_phi_deta2_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta2_norm_name,delta_phi_deta2_norm_p8_jetmet);
    if (delta_phi_deta2_norm_p8_jetmet == 0) { cout << delta_phi_deta2_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 
    
    TH1D *delta_phi_deta2_norm;
    delta_phi_deta2_norm = new TH1D("ak5PF_delta_phi_deta2_norm","Delta_phi_deta2 Norm;#Delta#phi [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma^{2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta2 Norm"<<endl; }
    calc_nomalization(delta_phi_deta2_norm_p6_jetmettau, delta_phi_deta2_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta2_norm, delta_phi_deta2_norm_p6_jetmettau, delta_phi_deta2_norm_p6_jetmet, delta_phi_deta2_norm_p8_jetmettau, delta_phi_deta2_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta2_norm, delta_phi_deta2_norm_p6_jetmettau, delta_phi_deta2_norm_p6_jetmet, delta_phi_deta2_norm_p8_jetmettau, delta_phi_deta2_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_norm, delta_phi_deta2_norm_p6_jetmettau, delta_phi_deta2_norm_p6_jetmet, delta_phi_deta2_norm_p8_jetmettau, delta_phi_deta2_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2_norm, delta_phi_deta2_norm_p6_jetmettau, delta_phi_deta2_norm_p6_jetmet, delta_phi_deta2_norm_p8_jetmettau, delta_phi_deta2_norm_p8_jetmet, final_errors, 25, detail);


//merging the delta phi deta3 distributions    
    TH1D *delta_phi_deta3_p6_jetmettau  = 0;
    TH1D *delta_phi_deta3_p6_jetmet  = 0;
    TH1D *delta_phi_deta3_p8_jetmettau  = 0;
    TH1D *delta_phi_deta3_p8_jetmet  = 0;
    TString delta_phi_deta3_name = "output_true_delta_phi_deta3";

    data_p6_jetmettau->GetObject(delta_phi_deta3_name,delta_phi_deta3_p6_jetmettau);
    if (delta_phi_deta3_p6_jetmettau == 0) { cout << delta_phi_deta3_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta3_name,delta_phi_deta3_p6_jetmet);
    if (delta_phi_deta3_p6_jetmet == 0) { cout << delta_phi_deta3_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta3_name,delta_phi_deta3_p8_jetmettau);
    if (delta_phi_deta3_p8_jetmettau == 0) { cout << delta_phi_deta3_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta3_name,delta_phi_deta3_p8_jetmet);
    if (delta_phi_deta3_p8_jetmet == 0) { cout << delta_phi_deta3_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    TH1D *delta_phi_deta3;
    delta_phi_deta3 = new TH1D("ak5PF_delta_phi_deta3","Delta_phi_deta3;|#Delta#phi| [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta3 "<<endl; }
    calc_nomalization(delta_phi_deta3_p6_jetmettau, delta_phi_deta3_p6_jetmet, weigth, detail);
    merge(delta_phi_deta3, delta_phi_deta3_p6_jetmettau, delta_phi_deta3_p6_jetmet, delta_phi_deta3_p8_jetmettau, delta_phi_deta3_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta3, delta_phi_deta3_p6_jetmettau, delta_phi_deta3_p6_jetmet, delta_phi_deta3_p8_jetmettau, delta_phi_deta3_p8_jetmet, output_path_plots, "merging_delta_phi_deta3", "top_left", true, detail);
    plot_histogram(delta_phi_deta3, delta_phi_deta3_p6_jetmettau, delta_phi_deta3_p6_jetmet, delta_phi_deta3_p8_jetmettau, delta_phi_deta3_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3, delta_phi_deta3_p6_jetmettau, delta_phi_deta3_p6_jetmet, delta_phi_deta3_p8_jetmettau, delta_phi_deta3_p8_jetmet, final_errors, 7, detail);


//merging the delta phi deta3 norm distributions    
    TH1D *delta_phi_deta3_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta3_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta3_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta3_norm_p8_jetmet  = 0;
    TString delta_phi_deta3_norm_name = "output_true_delta_phi_deta3_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta3_norm_name,delta_phi_deta3_norm_p6_jetmettau);
    if (delta_phi_deta3_norm_p6_jetmettau == 0) { cout << delta_phi_deta3_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta3_norm_name,delta_phi_deta3_norm_p6_jetmet);
    if (delta_phi_deta3_norm_p6_jetmet == 0) { cout << delta_phi_deta3_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta3_norm_name,delta_phi_deta3_norm_p8_jetmettau);
    if (delta_phi_deta3_norm_p8_jetmettau == 0) { cout << delta_phi_deta3_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta3_norm_name,delta_phi_deta3_norm_p8_jetmet);
    if (delta_phi_deta3_norm_p8_jetmet == 0) { cout << delta_phi_deta3_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    TH1D *delta_phi_deta3_norm;
    delta_phi_deta3_norm = new TH1D("ak5PF_delta_phi_deta3_norm","Delta_phi_deta3;|#Delta#phi| [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta3 Norm"<<endl; }
    calc_nomalization(delta_phi_deta3_norm_p6_jetmettau, delta_phi_deta3_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta3_norm, delta_phi_deta3_norm_p6_jetmettau, delta_phi_deta3_norm_p6_jetmet, delta_phi_deta3_norm_p8_jetmettau, delta_phi_deta3_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta3_norm, delta_phi_deta3_norm_p6_jetmettau, delta_phi_deta3_norm_p6_jetmet, delta_phi_deta3_norm_p8_jetmettau, delta_phi_deta3_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_norm, delta_phi_deta3_norm_p6_jetmettau, delta_phi_deta3_norm_p6_jetmet, delta_phi_deta3_norm_p8_jetmettau, delta_phi_deta3_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3_norm, delta_phi_deta3_norm_p6_jetmettau, delta_phi_deta3_norm_p6_jetmet, delta_phi_deta3_norm_p8_jetmettau, delta_phi_deta3_norm_p8_jetmet, final_errors, 26, detail);


//merging the delta phi deta4 distributions         
    TH1D *delta_phi_deta4_p6_jetmettau  = 0;
    TH1D *delta_phi_deta4_p6_jetmet  = 0;
    TH1D *delta_phi_deta4_p8_jetmettau  = 0;
    TH1D *delta_phi_deta4_p8_jetmet  = 0;
    TString delta_phi_deta4_name = "output_true_delta_phi_deta4";

    data_p6_jetmettau->GetObject(delta_phi_deta4_name,delta_phi_deta4_p6_jetmettau);
    if (delta_phi_deta4_p6_jetmettau == 0) { cout << delta_phi_deta4_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta4_name,delta_phi_deta4_p6_jetmet);
    if (delta_phi_deta4_p6_jetmet == 0) { cout << delta_phi_deta4_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta4_name,delta_phi_deta4_p8_jetmettau);
    if (delta_phi_deta4_p8_jetmettau == 0) { cout << delta_phi_deta4_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta4_name,delta_phi_deta4_p8_jetmet);
    if (delta_phi_deta4_p8_jetmet == 0) { cout << delta_phi_deta4_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    TH1D *delta_phi_deta4;
    delta_phi_deta4 = new TH1D("ak5PF_delta_phi_deta4","Delta_phi_deta4;|#Delta#phi| [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta4 "<<endl; }
    calc_nomalization(delta_phi_deta4_p6_jetmettau, delta_phi_deta4_p6_jetmet, weigth, detail);
    merge(delta_phi_deta4, delta_phi_deta4_p6_jetmettau, delta_phi_deta4_p6_jetmet, delta_phi_deta4_p8_jetmettau, delta_phi_deta4_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta4, delta_phi_deta4_p6_jetmettau, delta_phi_deta4_p6_jetmet, delta_phi_deta4_p8_jetmettau, delta_phi_deta4_p8_jetmet, output_path_plots, "merging_delta_phi_deta4", "top_left", true, detail);
    plot_histogram(delta_phi_deta4, delta_phi_deta4_p6_jetmettau, delta_phi_deta4_p6_jetmet, delta_phi_deta4_p8_jetmettau, delta_phi_deta4_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4, delta_phi_deta4_p6_jetmettau, delta_phi_deta4_p6_jetmet, delta_phi_deta4_p8_jetmettau, delta_phi_deta4_p8_jetmet, final_errors, 8, detail);


//merging the delta phi deta4 Norm distributions         
    TH1D *delta_phi_deta4_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta4_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta4_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta4_norm_p8_jetmet  = 0;
    TString delta_phi_deta4_norm_name = "output_true_delta_phi_deta4_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta4_norm_name,delta_phi_deta4_norm_p6_jetmettau);
    if (delta_phi_deta4_norm_p6_jetmettau == 0) { cout << delta_phi_deta4_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta4_norm_name,delta_phi_deta4_norm_p6_jetmet);
    if (delta_phi_deta4_norm_p6_jetmet == 0) { cout << delta_phi_deta4_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }    
    data_p8_jetmettau->GetObject(delta_phi_deta4_norm_name,delta_phi_deta4_norm_p8_jetmettau);
    if (delta_phi_deta4_norm_p8_jetmettau == 0) { cout << delta_phi_deta4_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta4_norm_name,delta_phi_deta4_norm_p8_jetmet);
    if (delta_phi_deta4_norm_p8_jetmet == 0) { cout << delta_phi_deta4_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; } 

    TH1D *delta_phi_deta4_norm;
    delta_phi_deta4_norm = new TH1D("ak5PF_delta_phi_deta4_norm","Delta_phi_deta4_norm;#Delta#phi [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma^{-2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta4 Norm"<<endl; }
    calc_nomalization(delta_phi_deta4_norm_p6_jetmettau, delta_phi_deta4_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta4_norm, delta_phi_deta4_norm_p6_jetmettau, delta_phi_deta4_norm_p6_jetmet, delta_phi_deta4_norm_p8_jetmettau, delta_phi_deta4_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta4_norm, delta_phi_deta4_norm_p6_jetmettau, delta_phi_deta4_norm_p6_jetmet, delta_phi_deta4_norm_p8_jetmettau, delta_phi_deta4_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_norm, delta_phi_deta4_norm_p6_jetmettau, delta_phi_deta4_norm_p6_jetmet, delta_phi_deta4_norm_p8_jetmettau, delta_phi_deta4_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4_norm, delta_phi_deta4_norm_p6_jetmettau, delta_phi_deta4_norm_p6_jetmet, delta_phi_deta4_norm_p8_jetmettau, delta_phi_deta4_norm_p8_jetmet, final_errors, 27, detail);



//merging the delta phi with gap deta1 distributions
    TH1D *delta_phi_deta1_gap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta1_gap_p6_jetmet  = 0;
    TH1D *delta_phi_deta1_gap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta1_gap_p8_jetmet  = 0;
    TString delta_phi_deta1_gap_name = "output_true_delta_phi_deta1_gap";

    data_p6_jetmettau->GetObject(delta_phi_deta1_gap_name,delta_phi_deta1_gap_p6_jetmettau);
    if (delta_phi_deta1_gap_p6_jetmettau == 0) { cout << delta_phi_deta1_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta1_gap_name,delta_phi_deta1_gap_p6_jetmet);
    if (delta_phi_deta1_gap_p6_jetmet == 0) { cout << delta_phi_deta1_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta1_gap_name,delta_phi_deta1_gap_p8_jetmettau);
    if (delta_phi_deta1_gap_p8_jetmettau == 0) { cout << delta_phi_deta1_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta1_gap_name,delta_phi_deta1_gap_p8_jetmet);
    if (delta_phi_deta1_gap_p8_jetmet == 0) { cout << delta_phi_deta1_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap = new TH1D("ak5PF_delta_phi_deta1_gap","Delta_phi_deta1_gap;|#Delta#phi| [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta1 gap"<<endl; }
    calc_nomalization(delta_phi_deta1_gap_p6_jetmettau, delta_phi_deta1_gap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta1_gap, delta_phi_deta1_gap_p6_jetmettau, delta_phi_deta1_gap_p6_jetmet, delta_phi_deta1_gap_p8_jetmettau, delta_phi_deta1_gap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta1_gap, delta_phi_deta1_gap_p6_jetmettau, delta_phi_deta1_gap_p6_jetmet, delta_phi_deta1_gap_p8_jetmettau, delta_phi_deta1_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_gap, delta_phi_deta1_gap_p6_jetmettau, delta_phi_deta1_gap_p6_jetmet, delta_phi_deta1_gap_p8_jetmettau, delta_phi_deta1_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1_gap, delta_phi_deta1_gap_p6_jetmettau, delta_phi_deta1_gap_p6_jetmet, delta_phi_deta1_gap_p8_jetmettau, delta_phi_deta1_gap_p8_jetmet, final_errors, 9, detail);


//merging the delta phi with gap deta1 Norm distributions
    TH1D *delta_phi_deta1_gap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta1_gap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta1_gap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta1_gap_norm_p8_jetmet  = 0;
    TString delta_phi_deta1_gap_norm_name = "output_true_delta_phi_deta1_gap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta1_gap_norm_name,delta_phi_deta1_gap_norm_p6_jetmettau);
    if (delta_phi_deta1_gap_norm_p6_jetmettau == 0) { cout << delta_phi_deta1_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta1_gap_norm_name,delta_phi_deta1_gap_norm_p6_jetmet);
    if (delta_phi_deta1_gap_norm_p6_jetmet == 0) { cout << delta_phi_deta1_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta1_gap_norm_name,delta_phi_deta1_gap_norm_p8_jetmettau);
    if (delta_phi_deta1_gap_norm_p8_jetmettau == 0) { cout << delta_phi_deta1_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta1_gap_norm_name,delta_phi_deta1_gap_norm_p8_jetmet);
    if (delta_phi_deta1_gap_norm_p8_jetmet == 0) { cout << delta_phi_deta1_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap_norm;
    delta_phi_deta1_gap_norm = new TH1D("ak5PF_delta_phi_deta1_gap_norm","Delta_phi_deta1_gap_norm;#Delta#phi [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma^{2}~#sigma{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta1 gap norm"<<endl; }
    calc_nomalization(delta_phi_deta1_gap_norm_p6_jetmettau, delta_phi_deta1_gap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_p6_jetmettau, delta_phi_deta1_gap_norm_p6_jetmet, delta_phi_deta1_gap_norm_p8_jetmettau, delta_phi_deta1_gap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_p6_jetmettau, delta_phi_deta1_gap_norm_p6_jetmet, delta_phi_deta1_gap_norm_p8_jetmettau, delta_phi_deta1_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_gap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_p6_jetmettau, delta_phi_deta1_gap_norm_p6_jetmet, delta_phi_deta1_gap_norm_p8_jetmettau, delta_phi_deta1_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_gap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_p6_jetmettau, delta_phi_deta1_gap_norm_p6_jetmet, delta_phi_deta1_gap_norm_p8_jetmettau, delta_phi_deta1_gap_norm_p8_jetmet, final_errors, 28, detail);


//merging the delta phi with gap deta2 distributions    
    TH1D *delta_phi_deta2_gap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta2_gap_p6_jetmet  = 0;
    TH1D *delta_phi_deta2_gap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta2_gap_p8_jetmet  = 0;
    TString delta_phi_deta2_gap_name = "output_true_delta_phi_deta2_gap";

    data_p6_jetmettau->GetObject(delta_phi_deta2_gap_name,delta_phi_deta2_gap_p6_jetmettau);
    if (delta_phi_deta2_gap_p6_jetmettau == 0) { cout << delta_phi_deta2_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta2_gap_name,delta_phi_deta2_gap_p6_jetmet);
    if (delta_phi_deta2_gap_p6_jetmet == 0) { cout << delta_phi_deta2_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta2_gap_name,delta_phi_deta2_gap_p8_jetmettau);
    if (delta_phi_deta2_gap_p8_jetmettau == 0) { cout << delta_phi_deta2_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta2_gap_name,delta_phi_deta2_gap_p8_jetmet);
    if (delta_phi_deta2_gap_p8_jetmet == 0) { cout << delta_phi_deta2_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap = new TH1D("ak5PF_delta_phi_deta2_gap","lala;|#Delta#phi| [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta2"<<endl; }
    calc_nomalization(delta_phi_deta2_gap_p6_jetmettau, delta_phi_deta2_gap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta2_gap, delta_phi_deta2_gap_p6_jetmettau, delta_phi_deta2_gap_p6_jetmet, delta_phi_deta2_gap_p8_jetmettau, delta_phi_deta2_gap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta2_gap, delta_phi_deta2_gap_p6_jetmettau, delta_phi_deta2_gap_p6_jetmet, delta_phi_deta2_gap_p8_jetmettau, delta_phi_deta2_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_gap, delta_phi_deta2_gap_p6_jetmettau, delta_phi_deta2_gap_p6_jetmet, delta_phi_deta2_gap_p8_jetmettau, delta_phi_deta2_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2_gap, delta_phi_deta2_gap_p6_jetmettau, delta_phi_deta2_gap_p6_jetmet, delta_phi_deta2_gap_p8_jetmettau, delta_phi_deta2_gap_p8_jetmet, final_errors, 10, detail);


//merging the delta phi with gap deta2 norm distributions    
    TH1D *delta_phi_deta2_gap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta2_gap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta2_gap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta2_gap_norm_p8_jetmet  = 0;
    TString delta_phi_deta2_gap_norm_name = "output_true_delta_phi_deta2_gap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta2_gap_norm_name,delta_phi_deta2_gap_norm_p6_jetmettau);
    if (delta_phi_deta2_gap_norm_p6_jetmettau == 0) { cout << delta_phi_deta2_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta2_gap_norm_name,delta_phi_deta2_gap_norm_p6_jetmet);
    if (delta_phi_deta2_gap_norm_p6_jetmet == 0) { cout << delta_phi_deta2_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta2_gap_norm_name,delta_phi_deta2_gap_norm_p8_jetmettau);
    if (delta_phi_deta2_gap_norm_p8_jetmettau == 0) { cout << delta_phi_deta2_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta2_gap_norm_name,delta_phi_deta2_gap_norm_p8_jetmet);
    if (delta_phi_deta2_gap_norm_p8_jetmet == 0) { cout << delta_phi_deta2_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap_norm;
    delta_phi_deta2_gap_norm = new TH1D("ak5PF_delta_phi_deta2_gap_norm","Delta Phi Deta2 Gap;#Delta#phi [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma^{2}~#sigma{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta2 Norm"<<endl; }
    calc_nomalization(delta_phi_deta2_gap_norm_p6_jetmettau, delta_phi_deta2_gap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_p6_jetmettau, delta_phi_deta2_gap_norm_p6_jetmet, delta_phi_deta2_gap_norm_p8_jetmettau, delta_phi_deta2_gap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_p6_jetmettau, delta_phi_deta2_gap_norm_p6_jetmet, delta_phi_deta2_gap_norm_p8_jetmettau, delta_phi_deta2_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_gap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_p6_jetmettau, delta_phi_deta2_gap_norm_p6_jetmet, delta_phi_deta2_gap_norm_p8_jetmettau, delta_phi_deta2_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_gap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_p6_jetmettau, delta_phi_deta2_gap_norm_p6_jetmet, delta_phi_deta2_gap_norm_p8_jetmettau, delta_phi_deta2_gap_norm_p8_jetmet, final_errors, 29, detail);


//merging the delta phi with gap deta3 distributions    
    TH1D *delta_phi_deta3_gap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta3_gap_p6_jetmet  = 0;
    TH1D *delta_phi_deta3_gap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta3_gap_p8_jetmet  = 0;
    TString delta_phi_deta3_gap_name = "output_true_delta_phi_deta3_gap";

    data_p6_jetmettau->GetObject(delta_phi_deta3_gap_name,delta_phi_deta3_gap_p6_jetmettau);
    if (delta_phi_deta3_gap_p6_jetmettau == 0) { cout << delta_phi_deta3_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta3_gap_name,delta_phi_deta3_gap_p6_jetmet);
    if (delta_phi_deta3_gap_p6_jetmet == 0) { cout << delta_phi_deta3_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta3_gap_name,delta_phi_deta3_gap_p8_jetmettau);
    if (delta_phi_deta3_gap_p8_jetmettau == 0) { cout << delta_phi_deta3_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta3_gap_name,delta_phi_deta3_gap_p8_jetmet);
    if (delta_phi_deta3_gap_p8_jetmet == 0) { cout << delta_phi_deta3_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap = new TH1D("ak5PF_delta_phi_deta3_gap","Delta_phi_deta3_gap;|#Delta#phi| [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta3"<<endl; }
    calc_nomalization(delta_phi_deta3_gap_p6_jetmettau, delta_phi_deta3_gap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta3_gap, delta_phi_deta3_gap_p6_jetmettau, delta_phi_deta3_gap_p6_jetmet, delta_phi_deta3_gap_p8_jetmettau, delta_phi_deta3_gap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta3_gap, delta_phi_deta3_gap_p6_jetmettau, delta_phi_deta3_gap_p6_jetmet, delta_phi_deta3_gap_p8_jetmettau, delta_phi_deta3_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_gap, delta_phi_deta3_gap_p6_jetmettau, delta_phi_deta3_gap_p6_jetmet, delta_phi_deta3_gap_p8_jetmettau, delta_phi_deta3_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3_gap, delta_phi_deta3_gap_p6_jetmettau, delta_phi_deta3_gap_p6_jetmet, delta_phi_deta3_gap_p8_jetmettau, delta_phi_deta3_gap_p8_jetmet, final_errors, 11, detail);


//merging the delta phi with gap deta3 norm distributions    
    TH1D *delta_phi_deta3_gap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta3_gap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta3_gap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta3_gap_norm_p8_jetmet  = 0;
    TString delta_phi_deta3_gap_norm_name = "output_true_delta_phi_deta3_gap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta3_gap_norm_name,delta_phi_deta3_gap_norm_p6_jetmettau);
    if (delta_phi_deta3_gap_norm_p6_jetmettau == 0) { cout << delta_phi_deta3_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta3_gap_norm_name,delta_phi_deta3_gap_norm_p6_jetmet);
    if (delta_phi_deta3_gap_norm_p6_jetmet == 0) { cout << delta_phi_deta3_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta3_gap_norm_name,delta_phi_deta3_gap_norm_p8_jetmettau);
    if (delta_phi_deta3_gap_norm_p8_jetmettau == 0) { cout << delta_phi_deta3_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta3_gap_norm_name,delta_phi_deta3_gap_norm_p8_jetmet);
    if (delta_phi_deta3_gap_norm_p8_jetmet == 0) { cout << delta_phi_deta3_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap_norm;
    delta_phi_deta3_gap_norm = new TH1D("ak5PF_delta_phi_deta3_gap_norm","Delta_phi_deta3_gap_norm;#Delta#phi [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma^{2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta3 Norm"<<endl; }
    calc_nomalization(delta_phi_deta3_gap_norm_p6_jetmettau, delta_phi_deta3_gap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_p6_jetmettau, delta_phi_deta3_gap_norm_p6_jetmet, delta_phi_deta3_gap_norm_p8_jetmettau, delta_phi_deta3_gap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_p6_jetmettau, delta_phi_deta3_gap_norm_p6_jetmet, delta_phi_deta3_gap_norm_p8_jetmettau, delta_phi_deta3_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_gap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_p6_jetmettau, delta_phi_deta3_gap_norm_p6_jetmet, delta_phi_deta3_gap_norm_p8_jetmettau, delta_phi_deta3_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_gap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_p6_jetmettau, delta_phi_deta3_gap_norm_p6_jetmet, delta_phi_deta3_gap_norm_p8_jetmettau, delta_phi_deta3_gap_norm_p8_jetmet, final_errors, 30, detail);


//merging the delta phi with gap deta4 distributions
    TH1D *delta_phi_deta4_gap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta4_gap_p6_jetmet  = 0;
    TH1D *delta_phi_deta4_gap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta4_gap_p8_jetmet  = 0;
    TString delta_phi_deta4_gap_name = "output_true_delta_phi_deta4_gap";

    data_p6_jetmettau->GetObject(delta_phi_deta4_gap_name,delta_phi_deta4_gap_p6_jetmettau);
    if (delta_phi_deta4_gap_p6_jetmettau == 0) { cout << delta_phi_deta4_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta4_gap_name,delta_phi_deta4_gap_p6_jetmet);
    if (delta_phi_deta4_gap_p6_jetmet == 0) { cout << delta_phi_deta4_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta4_gap_name,delta_phi_deta4_gap_p8_jetmettau);
    if (delta_phi_deta4_gap_p8_jetmettau == 0) { cout << delta_phi_deta4_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta4_gap_name,delta_phi_deta4_gap_p8_jetmet);
    if (delta_phi_deta4_gap_p8_jetmet == 0) { cout << delta_phi_deta4_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap = new TH1D("ak5PF_delta_phi_deta4_gap","Delta_phi_deta4_gap;|#Delta#phi| [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta4"<<endl; }
    calc_nomalization(delta_phi_deta4_gap_p6_jetmettau, delta_phi_deta4_gap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta4_gap, delta_phi_deta4_gap_p6_jetmettau, delta_phi_deta4_gap_p6_jetmet, delta_phi_deta4_gap_p8_jetmettau, delta_phi_deta4_gap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta4_gap, delta_phi_deta4_gap_p6_jetmettau, delta_phi_deta4_gap_p6_jetmet, delta_phi_deta4_gap_p8_jetmettau, delta_phi_deta4_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_gap, delta_phi_deta4_gap_p6_jetmettau, delta_phi_deta4_gap_p6_jetmet, delta_phi_deta4_gap_p8_jetmettau, delta_phi_deta4_gap_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4_gap, delta_phi_deta4_gap_p6_jetmettau, delta_phi_deta4_gap_p6_jetmet, delta_phi_deta4_gap_p8_jetmettau, delta_phi_deta4_gap_p8_jetmet, final_errors, 12, detail);


//merging the delta phi with gap deta4 norm distributions
    TH1D *delta_phi_deta4_gap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta4_gap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta4_gap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta4_gap_norm_p8_jetmet  = 0;
    TString delta_phi_deta4_gap_norm_name = "output_true_delta_phi_deta4_gap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta4_gap_norm_name,delta_phi_deta4_gap_norm_p6_jetmettau);
    if (delta_phi_deta4_gap_norm_p6_jetmettau == 0) { cout << delta_phi_deta4_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta4_gap_norm_name,delta_phi_deta4_gap_norm_p6_jetmet);
    if (delta_phi_deta4_gap_norm_p6_jetmet == 0) { cout << delta_phi_deta4_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta4_gap_norm_name,delta_phi_deta4_gap_norm_p8_jetmettau);
    if (delta_phi_deta4_gap_norm_p8_jetmettau == 0) { cout << delta_phi_deta4_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta4_gap_norm_name,delta_phi_deta4_gap_norm_p8_jetmet);
    if (delta_phi_deta4_gap_norm_p8_jetmet == 0) { cout << delta_phi_deta4_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap_norm;
    delta_phi_deta4_gap_norm = new TH1D("ak5PF_delta_phi_deta4_gap_norm","Delta_phi_deta4_gap_norm;#Delta#phi [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma^{2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta4 Norm"<<endl; }
    calc_nomalization(delta_phi_deta4_gap_norm_p6_jetmettau, delta_phi_deta4_gap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_p6_jetmettau, delta_phi_deta4_gap_norm_p6_jetmet, delta_phi_deta4_gap_norm_p8_jetmettau, delta_phi_deta4_gap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_p6_jetmettau, delta_phi_deta4_gap_norm_p6_jetmet, delta_phi_deta4_gap_norm_p8_jetmettau, delta_phi_deta4_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_gap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_p6_jetmettau, delta_phi_deta4_gap_norm_p6_jetmet, delta_phi_deta4_gap_norm_p8_jetmettau, delta_phi_deta4_gap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_gap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_p6_jetmettau, delta_phi_deta4_gap_norm_p6_jetmet, delta_phi_deta4_gap_norm_p8_jetmettau, delta_phi_deta4_gap_norm_p8_jetmet, final_errors, 31, detail);


//merging the delta phi without gap deta1 distributions    
    TH1D *delta_phi_deta1_nogap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta1_nogap_p6_jetmet  = 0;
    TH1D *delta_phi_deta1_nogap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta1_nogap_p8_jetmet  = 0;
    TString delta_phi_deta1_nogap_name = "output_true_delta_phi_deta1_nogap";

    data_p6_jetmettau->GetObject(delta_phi_deta1_nogap_name,delta_phi_deta1_nogap_p6_jetmettau);
    if (delta_phi_deta1_nogap_p6_jetmettau == 0) { cout << delta_phi_deta1_nogap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta1_nogap_name,delta_phi_deta1_nogap_p6_jetmet);
    if (delta_phi_deta1_nogap_p6_jetmet == 0) { cout << delta_phi_deta1_nogap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta1_nogap_name,delta_phi_deta1_nogap_p8_jetmettau);
    if (delta_phi_deta1_nogap_p8_jetmettau == 0) { cout << delta_phi_deta1_nogap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta1_nogap_name,delta_phi_deta1_nogap_p8_jetmet);
    if (delta_phi_deta1_nogap_p8_jetmet == 0) { cout << delta_phi_deta1_nogap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap = new TH1D("ak5PF_delta_phi_deta1_nogap","lala;|#Delta#phi| [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta1"<<endl; }
    calc_nomalization(delta_phi_deta1_nogap_p6_jetmettau, delta_phi_deta1_nogap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta1_nogap, delta_phi_deta1_nogap_p6_jetmettau, delta_phi_deta1_nogap_p6_jetmet, delta_phi_deta1_nogap_p8_jetmettau, delta_phi_deta1_nogap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta1_nogap, delta_phi_deta1_nogap_p6_jetmettau, delta_phi_deta1_nogap_p6_jetmet, delta_phi_deta1_nogap_p8_jetmettau, delta_phi_deta1_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_nogap, delta_phi_deta1_nogap_p6_jetmettau, delta_phi_deta1_nogap_p6_jetmet, delta_phi_deta1_nogap_p8_jetmettau, delta_phi_deta1_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1_nogap, delta_phi_deta1_nogap_p6_jetmettau, delta_phi_deta1_nogap_p6_jetmet, delta_phi_deta1_nogap_p8_jetmettau, delta_phi_deta1_nogap_p8_jetmet, final_errors, 13, detail);


//merging the delta phi without gap deta1 norm distributions    
    TH1D *delta_phi_deta1_nogap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta1_nogap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta1_nogap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta1_nogap_norm_p8_jetmet  = 0;
    TString delta_phi_deta1_nogap_norm_name = "output_true_delta_phi_deta1_nogap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta1_nogap_norm_name,delta_phi_deta1_nogap_norm_p6_jetmettau);
    if (delta_phi_deta1_nogap_norm_p6_jetmettau == 0) { cout << delta_phi_deta1_nogap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta1_nogap_norm_name,delta_phi_deta1_nogap_norm_p6_jetmet);
    if (delta_phi_deta1_nogap_norm_p6_jetmet == 0) { cout << delta_phi_deta1_nogap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta1_nogap_norm_name,delta_phi_deta1_nogap_norm_p8_jetmettau);
    if (delta_phi_deta1_nogap_norm_p8_jetmettau == 0) { cout << delta_phi_deta1_nogap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta1_nogap_norm_name,delta_phi_deta1_nogap_norm_p8_jetmet);
    if (delta_phi_deta1_nogap_norm_p8_jetmet == 0) { cout << delta_phi_deta1_nogap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap_norm;
    delta_phi_deta1_nogap_norm = new TH1D("ak5PF_delta_phi_deta1_nogap_norm","Delta Phi Deta1 Nogap Norm;#Delta#phi [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma^{2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta1 norm"<<endl; }
    calc_nomalization(delta_phi_deta1_nogap_norm_p6_jetmettau, delta_phi_deta1_nogap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_p6_jetmettau, delta_phi_deta1_nogap_norm_p6_jetmet, delta_phi_deta1_nogap_norm_p8_jetmettau, delta_phi_deta1_nogap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_p6_jetmettau, delta_phi_deta1_nogap_norm_p6_jetmet, delta_phi_deta1_nogap_norm_p8_jetmettau, delta_phi_deta1_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_nogap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_p6_jetmettau, delta_phi_deta1_nogap_norm_p6_jetmet, delta_phi_deta1_nogap_norm_p8_jetmettau, delta_phi_deta1_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta1_nogap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_p6_jetmettau, delta_phi_deta1_nogap_norm_p6_jetmet, delta_phi_deta1_nogap_norm_p8_jetmettau, delta_phi_deta1_nogap_norm_p8_jetmet, final_errors, 32, detail);


//merging the delta phi without gap deta2 distributions
    TH1D *delta_phi_deta2_nogap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta2_nogap_p6_jetmet  = 0;
    TH1D *delta_phi_deta2_nogap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta2_nogap_p8_jetmet  = 0;
    TString delta_phi_deta2_nogap_name = "output_true_delta_phi_deta2_nogap";

    data_p6_jetmettau->GetObject(delta_phi_deta2_nogap_name,delta_phi_deta2_nogap_p6_jetmettau);
    if (delta_phi_deta2_nogap_p6_jetmettau == 0) { cout << delta_phi_deta2_nogap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta2_nogap_name,delta_phi_deta2_nogap_p6_jetmet);
    if (delta_phi_deta2_nogap_p6_jetmet == 0) { cout << delta_phi_deta2_nogap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta2_nogap_name,delta_phi_deta2_nogap_p8_jetmettau);
    if (delta_phi_deta2_nogap_p8_jetmettau == 0) { cout << delta_phi_deta2_nogap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta2_nogap_name,delta_phi_deta2_nogap_p8_jetmet);
    if (delta_phi_deta2_nogap_p8_jetmet == 0) { cout << delta_phi_deta2_nogap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap = new TH1D("ak5PF_delta_phi_deta2_nogap","lala;|#Delta#phi| [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta2"<<endl; }
    calc_nomalization(delta_phi_deta2_nogap_p6_jetmettau, delta_phi_deta2_nogap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta2_nogap, delta_phi_deta2_nogap_p6_jetmettau, delta_phi_deta2_nogap_p6_jetmet, delta_phi_deta2_nogap_p8_jetmettau, delta_phi_deta2_nogap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta2_nogap, delta_phi_deta2_nogap_p6_jetmettau, delta_phi_deta2_nogap_p6_jetmet, delta_phi_deta2_nogap_p8_jetmettau, delta_phi_deta2_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_nogap, delta_phi_deta2_nogap_p6_jetmettau, delta_phi_deta2_nogap_p6_jetmet, delta_phi_deta2_nogap_p8_jetmettau, delta_phi_deta2_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2_nogap, delta_phi_deta2_nogap_p6_jetmettau, delta_phi_deta2_nogap_p6_jetmet, delta_phi_deta2_nogap_p8_jetmettau, delta_phi_deta2_nogap_p8_jetmet, final_errors, 14, detail);


//merging the delta phi without gap deta2 norm distributions
    TH1D *delta_phi_deta2_nogap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta2_nogap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta2_nogap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta2_nogap_norm_p8_jetmet  = 0;
    TString delta_phi_deta2_nogap_norm_name = "output_true_delta_phi_deta2_nogap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta2_nogap_norm_name,delta_phi_deta2_nogap_norm_p6_jetmettau);
    if (delta_phi_deta2_nogap_norm_p6_jetmettau == 0) { cout << delta_phi_deta2_nogap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta2_nogap_norm_name,delta_phi_deta2_nogap_norm_p6_jetmet);
    if (delta_phi_deta2_nogap_norm_p6_jetmet == 0) { cout << delta_phi_deta2_nogap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta2_nogap_norm_name,delta_phi_deta2_nogap_norm_p8_jetmettau);
    if (delta_phi_deta2_nogap_norm_p8_jetmettau == 0) { cout << delta_phi_deta2_nogap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta2_nogap_norm_name,delta_phi_deta2_nogap_norm_p8_jetmet);
    if (delta_phi_deta2_nogap_norm_p8_jetmet == 0) { cout << delta_phi_deta2_nogap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap_norm;
    delta_phi_deta2_nogap_norm = new TH1D("ak5PF_delta_phi_deta2_nogap_norm","Delta Phi Deta2 Nogap;#Delta#phi [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma^{2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta2 Norm"<<endl; }
    calc_nomalization(delta_phi_deta2_nogap_norm_p6_jetmettau, delta_phi_deta2_nogap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_p6_jetmettau, delta_phi_deta2_nogap_norm_p6_jetmet, delta_phi_deta2_nogap_norm_p8_jetmettau, delta_phi_deta2_nogap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_p6_jetmettau, delta_phi_deta2_nogap_norm_p6_jetmet, delta_phi_deta2_nogap_norm_p8_jetmettau, delta_phi_deta2_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_nogap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_p6_jetmettau, delta_phi_deta2_nogap_norm_p6_jetmet, delta_phi_deta2_nogap_norm_p8_jetmettau, delta_phi_deta2_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta2_nogap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_p6_jetmettau, delta_phi_deta2_nogap_norm_p6_jetmet, delta_phi_deta2_nogap_norm_p8_jetmettau, delta_phi_deta2_nogap_norm_p8_jetmet, final_errors, 33, detail);


//merging the delta phi without gap deta3 distributions
    TH1D *delta_phi_deta3_nogap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta3_nogap_p6_jetmet  = 0;
    TH1D *delta_phi_deta3_nogap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta3_nogap_p8_jetmet  = 0;
    TString delta_phi_deta3_nogap_name = "output_true_delta_phi_deta3_nogap";

    data_p6_jetmettau->GetObject(delta_phi_deta3_nogap_name,delta_phi_deta3_nogap_p6_jetmettau);
    if (delta_phi_deta3_nogap_p6_jetmettau == 0) { cout << delta_phi_deta3_nogap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta3_nogap_name,delta_phi_deta3_nogap_p6_jetmet);
    if (delta_phi_deta3_nogap_p6_jetmet == 0) { cout << delta_phi_deta3_nogap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta3_nogap_name,delta_phi_deta3_nogap_p8_jetmettau);
    if (delta_phi_deta3_nogap_p8_jetmettau == 0) { cout << delta_phi_deta3_nogap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta3_nogap_name,delta_phi_deta3_nogap_p8_jetmet);
    if (delta_phi_deta3_nogap_p8_jetmet == 0) { cout << delta_phi_deta3_nogap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap = new TH1D("ak5PF_delta_phi_deta3_nogap","lala;|#Delta#phi| [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta3"<<endl; }
    calc_nomalization(delta_phi_deta3_nogap_p6_jetmettau, delta_phi_deta3_nogap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta3_nogap, delta_phi_deta3_nogap_p6_jetmettau, delta_phi_deta3_nogap_p6_jetmet, delta_phi_deta3_nogap_p8_jetmettau, delta_phi_deta3_nogap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta3_nogap, delta_phi_deta3_nogap_p6_jetmettau, delta_phi_deta3_nogap_p6_jetmet, delta_phi_deta3_nogap_p8_jetmettau, delta_phi_deta3_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_nogap, delta_phi_deta3_nogap_p6_jetmettau, delta_phi_deta3_nogap_p6_jetmet, delta_phi_deta3_nogap_p8_jetmettau, delta_phi_deta3_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3_nogap, delta_phi_deta3_nogap_p6_jetmettau, delta_phi_deta3_nogap_p6_jetmet, delta_phi_deta3_nogap_p8_jetmettau, delta_phi_deta3_nogap_p8_jetmet, final_errors, 15, detail);


//merging the delta phi without gap deta3 norm distributions
    TH1D *delta_phi_deta3_nogap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta3_nogap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta3_nogap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta3_nogap_norm_p8_jetmet  = 0;
    TString delta_phi_deta3_nogap_norm_name = "output_true_delta_phi_deta3_nogap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta3_nogap_norm_name,delta_phi_deta3_nogap_norm_p6_jetmettau);
    if (delta_phi_deta3_nogap_norm_p6_jetmettau == 0) { cout << delta_phi_deta3_nogap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta3_nogap_norm_name,delta_phi_deta3_nogap_norm_p6_jetmet);
    if (delta_phi_deta3_nogap_norm_p6_jetmet == 0) { cout << delta_phi_deta3_nogap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta3_nogap_norm_name,delta_phi_deta3_nogap_norm_p8_jetmettau);
    if (delta_phi_deta3_nogap_norm_p8_jetmettau == 0) { cout << delta_phi_deta3_nogap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta3_nogap_norm_name,delta_phi_deta3_nogap_norm_p8_jetmet);
    if (delta_phi_deta3_nogap_norm_p8_jetmet == 0) { cout << delta_phi_deta3_nogap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap_norm;
    delta_phi_deta3_nogap_norm = new TH1D("ak5PF_delta_phi_deta3_nogap_norm","Delta Phi Deta3 Nogap Norm;#Delta#phi [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma^{2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta3 Norm"<<endl; }
    calc_nomalization(delta_phi_deta3_nogap_norm_p6_jetmettau, delta_phi_deta3_nogap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_p6_jetmettau, delta_phi_deta3_nogap_norm_p6_jetmet, delta_phi_deta3_nogap_norm_p8_jetmettau, delta_phi_deta3_nogap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_p6_jetmettau, delta_phi_deta3_nogap_norm_p6_jetmet, delta_phi_deta3_nogap_norm_p8_jetmettau, delta_phi_deta3_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_nogap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_p6_jetmettau, delta_phi_deta3_nogap_norm_p6_jetmet, delta_phi_deta3_nogap_norm_p8_jetmettau, delta_phi_deta3_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta3_nogap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_p6_jetmettau, delta_phi_deta3_nogap_norm_p6_jetmet, delta_phi_deta3_nogap_norm_p8_jetmettau, delta_phi_deta3_nogap_norm_p8_jetmet, final_errors, 34, detail);


//merging the delta phi without gap deta4 distributions
    TH1D *delta_phi_deta4_nogap_p6_jetmettau  = 0;
    TH1D *delta_phi_deta4_nogap_p6_jetmet  = 0;
    TH1D *delta_phi_deta4_nogap_p8_jetmettau  = 0;
    TH1D *delta_phi_deta4_nogap_p8_jetmet  = 0;
    TString delta_phi_deta4_nogap_name = "output_true_delta_phi_deta4_nogap";

    data_p6_jetmettau->GetObject(delta_phi_deta4_nogap_name,delta_phi_deta4_nogap_p6_jetmettau);
    if (delta_phi_deta4_nogap_p6_jetmettau == 0) { cout << delta_phi_deta4_nogap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta4_nogap_name,delta_phi_deta4_nogap_p6_jetmet);
    if (delta_phi_deta4_nogap_p6_jetmet == 0) { cout << delta_phi_deta4_nogap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta4_nogap_name,delta_phi_deta4_nogap_p8_jetmettau);
    if (delta_phi_deta4_nogap_p8_jetmettau == 0) { cout << delta_phi_deta4_nogap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta4_nogap_name,delta_phi_deta4_nogap_p8_jetmet);
    if (delta_phi_deta4_nogap_p8_jetmet == 0) { cout << delta_phi_deta4_nogap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap = new TH1D("ak5PF_delta_phi_deta4_nogap","lala;|#Delta#phi| [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta4"<<endl; }
    calc_nomalization(delta_phi_deta4_nogap_p6_jetmettau, delta_phi_deta4_nogap_p6_jetmet, weigth, detail);
    merge(delta_phi_deta4_nogap, delta_phi_deta4_nogap_p6_jetmettau, delta_phi_deta4_nogap_p6_jetmet, delta_phi_deta4_nogap_p8_jetmettau, delta_phi_deta4_nogap_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta4_nogap, delta_phi_deta4_nogap_p6_jetmettau, delta_phi_deta4_nogap_p6_jetmet, delta_phi_deta4_nogap_p8_jetmettau, delta_phi_deta4_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_nogap, delta_phi_deta4_nogap_p6_jetmettau, delta_phi_deta4_nogap_p6_jetmet, delta_phi_deta4_nogap_p8_jetmettau, delta_phi_deta4_nogap_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4_nogap, delta_phi_deta4_nogap_p6_jetmettau, delta_phi_deta4_nogap_p6_jetmet, delta_phi_deta4_nogap_p8_jetmettau, delta_phi_deta4_nogap_p8_jetmet, final_errors, 16, detail);


//merging the delta phi without gap deta4 norm distributions
    TH1D *delta_phi_deta4_nogap_norm_p6_jetmettau  = 0;
    TH1D *delta_phi_deta4_nogap_norm_p6_jetmet  = 0;
    TH1D *delta_phi_deta4_nogap_norm_p8_jetmettau  = 0;
    TH1D *delta_phi_deta4_nogap_norm_p8_jetmet  = 0;
    TString delta_phi_deta4_nogap_norm_name = "output_true_delta_phi_deta4_nogap_norm";

    data_p6_jetmettau->GetObject(delta_phi_deta4_nogap_norm_name,delta_phi_deta4_nogap_norm_p6_jetmettau);
    if (delta_phi_deta4_nogap_norm_p6_jetmettau == 0) { cout << delta_phi_deta4_nogap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_phi_deta4_nogap_norm_name,delta_phi_deta4_nogap_norm_p6_jetmet);
    if (delta_phi_deta4_nogap_norm_p6_jetmet == 0) { cout << delta_phi_deta4_nogap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_phi_deta4_nogap_norm_name,delta_phi_deta4_nogap_norm_p8_jetmettau);
    if (delta_phi_deta4_nogap_norm_p8_jetmettau == 0) { cout << delta_phi_deta4_nogap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_phi_deta4_nogap_norm_name,delta_phi_deta4_nogap_norm_p8_jetmet);
    if (delta_phi_deta4_nogap_norm_p8_jetmet == 0) { cout << delta_phi_deta4_nogap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap_norm;
    delta_phi_deta4_nogap_norm = new TH1D("ak5PF_delta_phi_deta4_nogap_norm","Delta Phi Deta4 NoGap;#Delta#phi [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma^{2}~#sigma^{-2}}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta4 Norm"<<endl; }
    calc_nomalization(delta_phi_deta4_nogap_norm_p6_jetmettau, delta_phi_deta4_nogap_norm_p6_jetmet, weigth, detail);
    merge(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_p6_jetmettau, delta_phi_deta4_nogap_norm_p6_jetmet, delta_phi_deta4_nogap_norm_p8_jetmettau, delta_phi_deta4_nogap_norm_p8_jetmet, weigth);
    plot_histogram(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_p6_jetmettau, delta_phi_deta4_nogap_norm_p6_jetmet, delta_phi_deta4_nogap_norm_p8_jetmettau, delta_phi_deta4_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_nogap_norm", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_p6_jetmettau, delta_phi_deta4_nogap_norm_p6_jetmet, delta_phi_deta4_nogap_norm_p8_jetmettau, delta_phi_deta4_nogap_norm_p8_jetmet, output_path_plots, "merging_delta_phi_deta4_nogap_norm_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_p6_jetmettau, delta_phi_deta4_nogap_norm_p6_jetmet, delta_phi_deta4_nogap_norm_p8_jetmettau, delta_phi_deta4_nogap_norm_p8_jetmet, final_errors, 35, detail);


//merging the leading eta star inside gap   
    TH1D *leading_eta_star_inside_gap_p6_jetmettau  = 0;
    TH1D *leading_eta_star_inside_gap_p6_jetmet  = 0;
    TH1D *leading_eta_star_inside_gap_p8_jetmettau  = 0;
    TH1D *leading_eta_star_inside_gap_p8_jetmet  = 0;
    TString leading_eta_star_inside_gap_name = "output_true_leading_eta_star_inside_gap";

    data_p6_jetmettau->GetObject(leading_eta_star_inside_gap_name,leading_eta_star_inside_gap_p6_jetmettau);
    if (leading_eta_star_inside_gap_p6_jetmettau == 0) { cout << leading_eta_star_inside_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(leading_eta_star_inside_gap_name,leading_eta_star_inside_gap_p6_jetmet);
    if (leading_eta_star_inside_gap_p6_jetmet == 0) { cout << leading_eta_star_inside_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(leading_eta_star_inside_gap_name,leading_eta_star_inside_gap_p8_jetmettau);
    if (leading_eta_star_inside_gap_p8_jetmettau == 0) { cout << leading_eta_star_inside_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(leading_eta_star_inside_gap_name,leading_eta_star_inside_gap_p8_jetmet);
    if (leading_eta_star_inside_gap_p8_jetmet == 0) { cout << leading_eta_star_inside_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap = new TH1D("ak5PF_leading_eta_star_inside_gap","Leading_eta_star_inside_gap;#eta* ;#frac{d#sigma}{d#eta*}", etastar_nbins, etastar_bins);
    leading_eta_star_inside_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_leading Eta Star Gap"<<endl; }
    calc_nomalization(leading_eta_star_inside_gap_p6_jetmettau, leading_eta_star_inside_gap_p6_jetmet, weigth, detail);
    merge(leading_eta_star_inside_gap, leading_eta_star_inside_gap_p6_jetmettau, leading_eta_star_inside_gap_p6_jetmet, leading_eta_star_inside_gap_p8_jetmettau, leading_eta_star_inside_gap_p8_jetmet, weigth);
    plot_histogram(leading_eta_star_inside_gap, leading_eta_star_inside_gap_p6_jetmettau, leading_eta_star_inside_gap_p6_jetmet, leading_eta_star_inside_gap_p8_jetmettau, leading_eta_star_inside_gap_p8_jetmet, output_path_plots, "merging_leading_eta_star_inside_gap", "bottom_middle", true, detail);
    plot_histogram(leading_eta_star_inside_gap, leading_eta_star_inside_gap_p6_jetmettau, leading_eta_star_inside_gap_p6_jetmet, leading_eta_star_inside_gap_p8_jetmettau, leading_eta_star_inside_gap_p8_jetmet, output_path_plots, "merging_leading_eta_star_inside_gap_linear", "bottom_left", false, detail);
    calc_error(leading_eta_star_inside_gap, leading_eta_star_inside_gap_p6_jetmettau, leading_eta_star_inside_gap_p6_jetmet, leading_eta_star_inside_gap_p8_jetmettau, leading_eta_star_inside_gap_p8_jetmet, final_errors, 17, detail);


//merging the leading eta star inside gap Norm   
    TH1D *leading_eta_star_inside_gap_norm_p6_jetmettau  = 0;
    TH1D *leading_eta_star_inside_gap_norm_p6_jetmet  = 0;
    TH1D *leading_eta_star_inside_gap_norm_p8_jetmettau  = 0;
    TH1D *leading_eta_star_inside_gap_norm_p8_jetmet  = 0;
    TString leading_eta_star_inside_gap_norm_name = "output_true_leading_eta_star_inside_gap_norm";

    data_p6_jetmettau->GetObject(leading_eta_star_inside_gap_norm_name,leading_eta_star_inside_gap_norm_p6_jetmettau);
    if (leading_eta_star_inside_gap_norm_p6_jetmettau == 0) { cout << leading_eta_star_inside_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(leading_eta_star_inside_gap_norm_name,leading_eta_star_inside_gap_norm_p6_jetmet);
    if (leading_eta_star_inside_gap_norm_p6_jetmet == 0) { cout << leading_eta_star_inside_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(leading_eta_star_inside_gap_norm_name,leading_eta_star_inside_gap_norm_p8_jetmettau);
    if (leading_eta_star_inside_gap_norm_p8_jetmettau == 0) { cout << leading_eta_star_inside_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(leading_eta_star_inside_gap_norm_name,leading_eta_star_inside_gap_norm_p8_jetmet);
    if (leading_eta_star_inside_gap_norm_p8_jetmet == 0) { cout << leading_eta_star_inside_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap_norm;
    leading_eta_star_inside_gap_norm = new TH1D("ak5PF_leading_eta_star_inside_gap_norm","Leading_eta_star_inside_gap_norm;#eta* ;#frac{d#sigma~#sigma{-1}}{d#eta*}", etastar_nbins, etastar_bins);
    leading_eta_star_inside_gap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_leading Eta Star Gap Norm"<<endl; }
    calc_nomalization(leading_eta_star_inside_gap_norm_p6_jetmettau, leading_eta_star_inside_gap_norm_p6_jetmet, weigth, detail);
    merge(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_p6_jetmettau, leading_eta_star_inside_gap_norm_p6_jetmet, leading_eta_star_inside_gap_norm_p8_jetmettau, leading_eta_star_inside_gap_norm_p8_jetmet, weigth);
    plot_histogram(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_p6_jetmettau, leading_eta_star_inside_gap_norm_p6_jetmet, leading_eta_star_inside_gap_norm_p8_jetmettau, leading_eta_star_inside_gap_norm_p8_jetmet, output_path_plots, "merging_leading_eta_star_inside_gap_norm", "bottom_middle", true, detail);
    plot_histogram(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_p6_jetmettau, leading_eta_star_inside_gap_norm_p6_jetmet, leading_eta_star_inside_gap_norm_p8_jetmettau, leading_eta_star_inside_gap_norm_p8_jetmet, output_path_plots, "merging_leading_eta_star_inside_gap_norm_linear", "bottom_left", false, detail);
    calc_error(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_p6_jetmettau, leading_eta_star_inside_gap_norm_p6_jetmet, leading_eta_star_inside_gap_norm_p8_jetmettau, leading_eta_star_inside_gap_norm_p8_jetmet, final_errors, 36, detail);


//merging the delta eta outside gap distribution
    TH1D *delta_eta_outside_gap_p6_jetmettau  = 0;
    TH1D *delta_eta_outside_gap_p6_jetmet  = 0;
    TH1D *delta_eta_outside_gap_p8_jetmettau  = 0;
    TH1D *delta_eta_outside_gap_p8_jetmet  = 0;
    TString delta_eta_outside_gap_name = "output_true_delta_eta_outside_gap";

    data_p6_jetmettau->GetObject(delta_eta_outside_gap_name,delta_eta_outside_gap_p6_jetmettau);
    if (delta_eta_outside_gap_p6_jetmettau == 0) { cout << delta_eta_outside_gap_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_eta_outside_gap_name,delta_eta_outside_gap_p6_jetmet);
    if (delta_eta_outside_gap_p6_jetmet == 0) { cout << delta_eta_outside_gap_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_eta_outside_gap_name,delta_eta_outside_gap_p8_jetmettau);
    if (delta_eta_outside_gap_p8_jetmettau == 0) { cout << delta_eta_outside_gap_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_eta_outside_gap_name,delta_eta_outside_gap_p8_jetmet);
    if (delta_eta_outside_gap_p8_jetmet == 0) { cout << delta_eta_outside_gap_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap = new TH1D("ak5PF_delta_eta_outside_gap","delta_eta_outside_gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta}", deta_out_nbins, deta_out_bins);
    delta_eta_outside_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta eta outside gap"<<endl; }
    calc_nomalization(delta_eta_outside_gap_p6_jetmettau, delta_eta_outside_gap_p6_jetmet, weigth, detail);
    merge(delta_eta_outside_gap, delta_eta_outside_gap_p6_jetmettau, delta_eta_outside_gap_p6_jetmet, delta_eta_outside_gap_p8_jetmettau, delta_eta_outside_gap_p8_jetmet, weigth);
    plot_histogram(delta_eta_outside_gap, delta_eta_outside_gap_p6_jetmettau, delta_eta_outside_gap_p6_jetmet, delta_eta_outside_gap_p8_jetmettau, delta_eta_outside_gap_p8_jetmet, output_path_plots, "merging_delta_eta_outside_gap", "bottom_left", true, detail);
    plot_histogram(delta_eta_outside_gap, delta_eta_outside_gap_p6_jetmettau, delta_eta_outside_gap_p6_jetmet, delta_eta_outside_gap_p8_jetmettau, delta_eta_outside_gap_p8_jetmet, output_path_plots, "merging_delta_eta_outside_gap_linear", "top_right", false, detail);
    calc_error(delta_eta_outside_gap, delta_eta_outside_gap_p6_jetmettau, delta_eta_outside_gap_p6_jetmet, delta_eta_outside_gap_p8_jetmettau, delta_eta_outside_gap_p8_jetmet, final_errors, 18, detail);


//merging the delta eta outside gap norm distribution
    TH1D *delta_eta_outside_gap_norm_p6_jetmettau  = 0;
    TH1D *delta_eta_outside_gap_norm_p6_jetmet  = 0;
    TH1D *delta_eta_outside_gap_norm_p8_jetmettau  = 0;
    TH1D *delta_eta_outside_gap_norm_p8_jetmet  = 0;
    TString delta_eta_outside_gap_norm_name = "output_true_delta_eta_outside_gap_norm";

    data_p6_jetmettau->GetObject(delta_eta_outside_gap_norm_name,delta_eta_outside_gap_norm_p6_jetmettau);
    if (delta_eta_outside_gap_norm_p6_jetmettau == 0) { cout << delta_eta_outside_gap_norm_name << " for JetMETTau_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p6_jetmet->GetObject(delta_eta_outside_gap_norm_name,delta_eta_outside_gap_norm_p6_jetmet);
    if (delta_eta_outside_gap_norm_p6_jetmet == 0) { cout << delta_eta_outside_gap_norm_name << " for JetMET_2010A Pythia 6 Z2* not found!" << endl; return; }
    data_p8_jetmettau->GetObject(delta_eta_outside_gap_norm_name,delta_eta_outside_gap_norm_p8_jetmettau);
    if (delta_eta_outside_gap_norm_p8_jetmettau == 0) { cout << delta_eta_outside_gap_norm_name << " for JetMETTau_2010A Pythia 8 4C not found!" << endl; return; }
    data_p8_jetmet->GetObject(delta_eta_outside_gap_norm_name,delta_eta_outside_gap_norm_p8_jetmet);
    if (delta_eta_outside_gap_norm_p8_jetmet == 0) { cout << delta_eta_outside_gap_norm_name << " for JetMET_2010A Pythia 8 4C not found!" << endl; return; }

    TH1D *delta_eta_outside_gap_norm;
    delta_eta_outside_gap_norm = new TH1D("ak5PF_delta_eta_outside_gap_norm","delta_eta_outside_gap_norm;#Delta#eta;#frac{d#sigma~#sigma^{-1}}{d#Delta#eta}", deta_out_nbins, deta_out_bins);
    delta_eta_outside_gap_norm->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta eta outside gap norm"<<endl; }
    calc_nomalization(delta_eta_outside_gap_norm_p6_jetmettau, delta_eta_outside_gap_norm_p6_jetmet, weigth, detail);
    merge(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_p6_jetmettau, delta_eta_outside_gap_norm_p6_jetmet, delta_eta_outside_gap_norm_p8_jetmettau, delta_eta_outside_gap_norm_p8_jetmet, weigth);
    plot_histogram(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_p6_jetmettau, delta_eta_outside_gap_norm_p6_jetmet, delta_eta_outside_gap_norm_p8_jetmettau, delta_eta_outside_gap_norm_p8_jetmet, output_path_plots, "merging_delta_eta_outside_gap_norm", "bottom_left", true, detail);
    plot_histogram(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_p6_jetmettau, delta_eta_outside_gap_norm_p6_jetmet, delta_eta_outside_gap_norm_p8_jetmettau, delta_eta_outside_gap_norm_p8_jetmet, output_path_plots, "merging_delta_eta_outside_gap_norm_linear", "top_right", false, detail);
    calc_error(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_p6_jetmettau, delta_eta_outside_gap_norm_p6_jetmet, delta_eta_outside_gap_norm_p8_jetmettau, delta_eta_outside_gap_norm_p8_jetmet, final_errors, 37, detail);


/*
//merging the delta eta distributions
    TH1D *delta_eta_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_eta");
    TH1D *delta_eta_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_eta");
    TH1D *delta_eta_jet = (TH1D*) data_jet->Get("ak5PF_delta_eta");

    TH1D *delta_eta;
    delta_eta = new TH1D("ak5PF_delta_eta","Delta_eta;#Delta#eta;#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Eta"<<endl; }
    calc_nomalization(delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, weigth, detail);
    merge(delta_eta, delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, weigth);
    plot_histogram(delta_eta, delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, output_path_plots, "merging_delta_eta", "bottom_left");
    calc_error(delta_eta, delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, final_errors, 23, detail);

//merging the delta eta with gap distributions
    TH1D *delta_eta_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_eta_gap");
    TH1D *delta_eta_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_eta_gap");
    TH1D *delta_eta_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_eta_gap");

    TH1D *delta_eta_gap;
    delta_eta_gap = new TH1D("ak5PF_delta_eta_gap","Delta_eta_gap;#Delta#eta^{jet-veto};#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Eta Gap"<<endl; }
    calc_nomalization(delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, weigth, detail);
    merge(delta_eta_gap, delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, weigth);
    plot_histogram(delta_eta_gap, delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, output_path_plots, "merging_delta_eta_gap", "bottom_left");
    calc_error(delta_eta_gap, delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, final_errors, 24, detail);

//merging the delta eta without gap distributions
    TH1D *delta_eta_nogap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_eta_nogap");
    TH1D *delta_eta_nogap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_eta_nogap");
    TH1D *delta_eta_nogap_jet = (TH1D*) data_jet->Get("ak5PF_delta_eta_nogap");

    TH1D *delta_eta_nogap;
    delta_eta_nogap = new TH1D("ak5PF_delta_eta_nogap","Delta_eta_nogap;#Delta#eta^{requiring a jet};#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Eta nogap"<<endl; }
    calc_nomalization(delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, weigth, detail);
    merge(delta_eta_nogap, delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, weigth);
    plot_histogram(delta_eta_nogap, delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, output_path_plots, "merging_delta_eta_nogap", "bottom_left");
    calc_error(delta_eta_nogap, delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, final_errors, 25, detail);

//merging the delta eta without gap distributions
    TH1D *vertex_selected_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_vertex_selected");
    TH1D *vertex_selected_jetmet = (TH1D*) data_jetmet->Get("ak5PF_vertex_selected");
    TH1D *vertex_selected_jet = (TH1D*) data_jet->Get("ak5PF_vertex_selected");

    TH1D *vertex_selected;
    vertex_selected = new TH1D("ak5PF_vertex_selected","Vertex multiplicity in selected events;Multiplicity;#frac{d#sigma}{dN} [pb]", 15, 0, 15);
    vertex_selected->Sumw2();
    
    if (detail) { cout<<"ak5PF_vertex_selected"<<endl; }
    calc_nomalization(vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, weigth, detail);
    merge(vertex_selected, vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, weigth);
    plot_histogram(vertex_selected, vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, output_path_plots, "merging_vertex_selected", "bottom_left", detail);
    calc_error(vertex_selected, vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, final_errors, 26, detail);
*/


//output the error variation
    if (detail) { cout<<"Display the final errors..."<<endl; }
    if (disp_errors) { show_error_variation2(final_errors); }

//creating the output file
    TFile data_output( output_rootfile.c_str() , "RECREATE");

//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file "<<output_rootfile<<" ..."<<endl; }
    delta_phi->Write();
    delta_phi_deta1->Write();
    delta_phi_deta2->Write();
    delta_phi_deta3->Write();
    delta_phi_deta4->Write();
    delta_phi_gap->Write();
    delta_phi_deta1_gap->Write();
    delta_phi_deta2_gap->Write();
    delta_phi_deta3_gap->Write();
    delta_phi_deta4_gap->Write();
    delta_phi_nogap->Write();
    delta_phi_deta1_nogap->Write();
    delta_phi_deta2_nogap->Write();
    delta_phi_deta3_nogap->Write();
    delta_phi_deta4_nogap->Write();
    leading_pt_inside_gap->Write();
    leading_pt_outside_gap->Write();
    leading_eta_star_inside_gap->Write();
    delta_eta_outside_gap->Write();

    delta_phi_norm->Write();
    delta_phi_deta1_norm->Write();
    delta_phi_deta2_norm->Write();
    delta_phi_deta3_norm->Write();
    delta_phi_deta4_norm->Write();
    delta_phi_gap_norm->Write();
    delta_phi_deta1_gap_norm->Write();
    delta_phi_deta2_gap_norm->Write();
    delta_phi_deta3_gap_norm->Write();
    delta_phi_deta4_gap_norm->Write();
    delta_phi_nogap_norm->Write();
    delta_phi_deta1_nogap_norm->Write();
    delta_phi_deta2_nogap_norm->Write();
    delta_phi_deta3_nogap_norm->Write();
    delta_phi_deta4_nogap_norm->Write();
    leading_pt_inside_gap_norm->Write();
    leading_pt_outside_gap_norm->Write();
    leading_eta_star_inside_gap_norm->Write();
    delta_eta_outside_gap_norm->Write();
    if (detail) { cout<<"Writing on "<<output_rootfile<<" was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data_p6_jetmettau->Close();
    data_p6_jetmet->Close();
    data_p8_jetmettau->Close();
    data_p8_jetmet->Close();
    data_output.Close();

    if (detail) { cout<<"Done!"<<endl; }
       
}
