// Pedro Cipriano, Nov 2013
// DESY, CMS
// Last Update: 08 Nov 2013
//
// compute_uncorr_uncertainty()

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TPad.h>
#include <TString.h>
#include <TF1.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "common_methods.h"

using namespace std;


void calc_uncorr_uncertainty(TH1 *uncorr, TH1 *stat, TH1 *model, TH1 *jer, double *uncorr_unc, int index = 0, string path = "../output/" , string fileout = "test", string legend_position = "top_left", bool detail = false)
{
//initialize variables
double max = 0.0;
double min = 0.0;
double tot = 0.0;
double ave = 0.0;

//set base values for the uncertainties
double val_stat, val_model, val_jer;
double val_pileup = 0.01;
double val_uncorr = 0.0;

//clonning the histograms
TH1D *pileup =  (TH1D*) stat->Clone();
TH1D *stat_plot =  (TH1D*) stat->Clone();

//loop over the bins
for (int i = 1; i <= uncorr->GetNbinsX();i++)
	{
	val_stat = stat->GetBinError(i)/stat->GetBinContent(i);
	val_model = model->GetBinContent(i);
        val_jer = jer->GetBinContent(i);
	val_uncorr = sqrt(val_stat * val_stat +  val_model * val_model + val_pileup * val_pileup + val_jer * val_jer);
	uncorr->SetBinContent(i,val_uncorr);
	pileup->SetBinContent(i,val_pileup);
	stat_plot->SetBinContent(i,val_stat);
        if (val_uncorr > max) {max = val_uncorr;}
        tot = tot + val_uncorr;
        if (val_uncorr < min || i == 1) { min = val_uncorr;}
        if (detail) { cout<< "Bin " << i << " Stat = " << val_stat << ", Model = " << val_model << ", Pileup = " << val_pileup << ", JER = " << val_jer <<" and Uncorrelated total = " << val_uncorr << endl; }
	}

//calculates the average model uncertainty
    ave = tot/uncorr->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    uncorr_unc[index*3+0] = ave*100;
    uncorr_unc[index*3+1] = min*100;
    uncorr_unc[index*3+2] = max*100;

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
// declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

//calculate the plooting range
    double min_plot = 0.0;
    double max_plot = 0.55;

//plooting
    uncorr->SetMaximum(max_plot);
    uncorr->SetMinimum(min_plot);
    format_histogram(uncorr, 1, 1);
    uncorr->Draw("hist");
    format_histogram(stat_plot, 2, 2);
    stat_plot->Draw("hist same");
    format_histogram(model, 4, 3);
    model->Draw("hist same");
    format_histogram(pileup, 6, 4);
    pileup->Draw("hist same");
    format_histogram(jer, 7, 5);
    jer->Draw("hist same");


//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 5, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(uncorr,"Unorrelated Uncertainty","l");
    leg01->AddEntry(stat_plot,"Statistical Uncertainty","l");
    leg01->AddEntry(model,"Model Depence Uncertainty","l");
    leg01->AddEntry(pileup,"Pileup Uncertainty","l");
    leg01->AddEntry(jer,"Jet Energy Resolution Uncertainty","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

//print the plot
    print_plots(c01, path, fileout);

}


void show_uncorr_uncertainties(double *uncorr_unc)
{
//shows the computed jes uncertainties
    cout<<" "<<endl;
    cout<<"Total Uncertainty"<<endl;
    cout<<"Observable                  Average  Minimum  Maximum"<<endl;
    cout<<"Delta phi                   "<<uncorr_unc[0]<<"  "<<uncorr_unc[1]<<" "<<uncorr_unc[2]<<endl;
    cout<<"Delta phi deta1             "<<uncorr_unc[3]<<"  "<<uncorr_unc[4]<<"  "<<uncorr_unc[5]<<endl;
    cout<<"Delta phi deta2             "<<uncorr_unc[6]<<"   "<<uncorr_unc[7]<<" "<<uncorr_unc[8]<<endl;
    cout<<"Delta phi deta3             "<<uncorr_unc[9]<<"  "<<uncorr_unc[10]<<" "<<uncorr_unc[11]<<endl;
    cout<<"Delta phi deta4             "<<uncorr_unc[12]<<"   "<<uncorr_unc[13]<<" "<<uncorr_unc[14]<<endl;
    cout<<"Delta phi gap               "<<uncorr_unc[15]<<"  "<<uncorr_unc[16]<<" "<<uncorr_unc[17]<<endl;
    cout<<"Delta phi deta1 gap         "<<uncorr_unc[18]<<"  "<<uncorr_unc[19]<<"  "<<uncorr_unc[20]<<endl;
    cout<<"Delta phi deta2 gap         "<<uncorr_unc[21]<<"  "<<uncorr_unc[22]<<" "<<uncorr_unc[23]<<endl;
    cout<<"Delta phi deta3 gap         "<<uncorr_unc[24]<<"  "<<uncorr_unc[25]<<" "<<uncorr_unc[26]<<endl;
    cout<<"Delta phi deta4 gap         "<<uncorr_unc[27]<<"  "<<uncorr_unc[28]<<" "<<uncorr_unc[29]<<endl;
    cout<<"Delta phi nogap             "<<uncorr_unc[30]<<"  "<<uncorr_unc[31]<<" "<<uncorr_unc[32]<<endl;
    cout<<"Delta phi deta1 nogap       "<<uncorr_unc[33]<<"  "<<uncorr_unc[34]<<"  "<<uncorr_unc[35]<<endl;
    cout<<"Delta phi deta2 nogap       "<<uncorr_unc[36]<<"   "<<uncorr_unc[37]<<"  "<<uncorr_unc[38]<<endl;
    cout<<"Delta phi deta3 nogap       "<<uncorr_unc[39]<<"  "<<uncorr_unc[40]<<"  "<<uncorr_unc[41]<<endl;
    cout<<"Delta phi deta4 nogap       "<<uncorr_unc[42]<<"  "<<uncorr_unc[43]<<"  "<<uncorr_unc[44]<<endl;
    cout<<"Leading pT inside gap       "<<uncorr_unc[45]<<"  "<<uncorr_unc[46]<<" "<<uncorr_unc[47]<<endl;
    cout<<"Leading eta star inside gap "<<uncorr_unc[48]<<"  "<<uncorr_unc[49]<<"  "<<uncorr_unc[50]<<endl;
    cout<<"Delta eta outside gap       "<<uncorr_unc[51]<<"  "<<uncorr_unc[52]<<" "<<uncorr_unc[53]<<endl;
    cout<<"Leading pT outside gap      "<<uncorr_unc[54]<<"  "<<uncorr_unc[55]<<" "<<uncorr_unc[56]<<endl;
}


void uncorr_uncertainty(string path_stat, string path_model, string path_jer, string uncorr_uncertainty, string label_out, string output_path_plots, bool detail = false, bool disp_uncertainty = true, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Total Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for Stat:  " << path_stat << endl; }
    if (detail) { cout << "Input path for Model: " << path_model << endl; }
    if (detail) { cout << "Input path for JER:   " << path_jer << endl; }
    if (detail) { cout << "Label Out:            " << label_out << endl; }
    if (detail) { cout << "Output path:          " << uncorr_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:    " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:         " << detail << endl; }
    if (detail) { cout << "Display Results:      " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:            " << test << endl; }

//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_stat = new TFile( path_stat.c_str() );
    TFile *data_model = new TFile( path_model.c_str() );
    TFile *data_jer = new TFile( path_jer.c_str() );

//histogram bins
int in_nbins = 9;
int out_nbins = 9;

double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//declaring and initializing the uncertainties array
    double uncorr_unc[19*3];
    double uncorr_unc_norm[19*3];

    for (int i=0; i<= 19*3-1;i++)
    {
    uncorr_unc[i] = 0.0;
    uncorr_unc_norm[i] = 0.0;
    }


//compute uncorrelated uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_stat = 0;
    TH1D *delta_phi_model = 0;
    TH1D *delta_phi_jer = 0;
    TString delta_phi_name_out = label_out + "delta_phi";

    data_stat->GetObject("ak5PF_delta_phi",delta_phi_stat);
    if (delta_phi_stat == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi",delta_phi_model);
    if (delta_phi_model == 0) { cout << "merged_model_unc_delta_phi not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi",delta_phi_jer);
    if (delta_phi_jer == 0) { cout << "merged_jer_unc_delta_phi not found!" << endl; return; }
    
    TH1D *delta_phi;
    delta_phi =  new TH1D(delta_phi_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi, delta_phi_stat, delta_phi_model, delta_phi_jer, uncorr_unc, 0, output_path_plots, label_out + "delta_phi", "top_right", detail);


//compute uncorrelated uncertainty for delta phi norm distribution
    if (detail) { cout<<"Delta phi Norm"<<endl; }

    TH1D *delta_phi_norm_stat = 0;
    TH1D *delta_phi_norm_model = 0;
    TH1D *delta_phi_norm_jer = 0;
    TString delta_phi_norm_name_out = label_out + "delta_phi_norm";

    data_stat->GetObject("ak5PF_delta_phi_norm",delta_phi_norm_stat);
    if (delta_phi_norm_stat == 0) { cout << "ak5PF_delta_phi_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_norm",delta_phi_norm_model);
    if (delta_phi_norm_model == 0) { cout << "merged_model_unc_delta_phi_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_norm",delta_phi_norm_jer);
    if (delta_phi_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_norm not found!" << endl; return; }
    
    TH1D *delta_phi_norm;
    delta_phi_norm =  new TH1D(delta_phi_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_norm, delta_phi_norm_stat, delta_phi_norm_model, delta_phi_norm_jer, uncorr_unc_norm, 0, output_path_plots, label_out + "delta_phi_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi Deta1"<<endl; }

    TH1D *delta_phi_deta1_stat = 0;
    TH1D *delta_phi_deta1_model = 0;
    TH1D *delta_phi_deta1_jer = 0;
    TString delta_phi_deta1_name_out = label_out + "delta_phi_deta1";

    data_stat->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_stat);
    if (delta_phi_deta1_stat == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1",delta_phi_deta1_model);
    if (delta_phi_deta1_model == 0) { cout << "merged_model_unc_delta_phi_deta1 not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta1",delta_phi_deta1_jer);
    if (delta_phi_deta1_jer == 0) { cout << "merged_jer_unc_delta_phi_deta1 not found!" << endl; return; }
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D(delta_phi_deta1_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta1, delta_phi_deta1_stat, delta_phi_deta1_model, delta_phi_deta1_jer, uncorr_unc, 1, output_path_plots, label_out + "delta_phi_deta1", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta1 norm distribution
    if (detail) { cout<<"Delta phi Deta1 Norm"<<endl; }

    TH1D *delta_phi_deta1_norm_stat = 0;
    TH1D *delta_phi_deta1_norm_model = 0;
    TH1D *delta_phi_deta1_norm_jer = 0;
    TString delta_phi_deta1_norm_name_out = label_out + "delta_phi_deta1_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta1_norm",delta_phi_deta1_norm_stat);
    if (delta_phi_deta1_norm_stat == 0) { cout << "ak5PF_delta_phi_deta1_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_norm",delta_phi_deta1_norm_model);
    if (delta_phi_deta1_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta1_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta1_norm",delta_phi_deta1_norm_jer);
    if (delta_phi_deta1_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta1_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_norm;
    delta_phi_deta1_norm =  new TH1D(delta_phi_deta1_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta1_norm, delta_phi_deta1_norm_stat, delta_phi_deta1_norm_model, delta_phi_deta1_norm_jer, uncorr_unc_norm, 1, output_path_plots, label_out + "delta_phi_deta1_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi Deta2"<<endl; }

    TH1D *delta_phi_deta2_stat = 0;
    TH1D *delta_phi_deta2_model = 0;
    TH1D *delta_phi_deta2_jer = 0;
    TString delta_phi_deta2_name_out = label_out + "delta_phi_deta2";

    data_stat->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_stat);
    if (delta_phi_deta2_stat == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2",delta_phi_deta2_model);
    if (delta_phi_deta2_model == 0) { cout << "merged_model_unc_delta_phi_deta2 not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta2",delta_phi_deta2_jer);
    if (delta_phi_deta2_jer == 0) { cout << "merged_jer_unc_delta_phi_deta2 not found!" << endl; return; }
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D(delta_phi_deta2_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta2, delta_phi_deta2_stat, delta_phi_deta2_model, delta_phi_deta2_jer, uncorr_unc, 2, output_path_plots, label_out + "delta_phi_deta2", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta2 norm distribution
    if (detail) { cout<<"Delta phi Deta2 Norm"<<endl; }

    TH1D *delta_phi_deta2_norm_stat = 0;
    TH1D *delta_phi_deta2_norm_model = 0;
    TH1D *delta_phi_deta2_norm_jer = 0;
    TString delta_phi_deta2_norm_name_out = label_out + "delta_phi_deta2_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta2_norm",delta_phi_deta2_norm_stat);
    if (delta_phi_deta2_norm_stat == 0) { cout << "ak5PF_delta_phi_deta2_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_norm",delta_phi_deta2_norm_model);
    if (delta_phi_deta2_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta2_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta2_norm",delta_phi_deta2_norm_jer);
    if (delta_phi_deta2_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta2_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_norm;
    delta_phi_deta2_norm =  new TH1D(delta_phi_deta2_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta2_norm, delta_phi_deta2_norm_stat, delta_phi_deta2_norm_model, delta_phi_deta2_norm_jer, uncorr_unc_norm, 2, output_path_plots, label_out + "delta_phi_deta2_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi Deta3"<<endl; }

    TH1D *delta_phi_deta3_stat = 0;
    TH1D *delta_phi_deta3_model = 0;
    TH1D *delta_phi_deta3_jer = 0;
    TString delta_phi_deta3_name_out = label_out + "delta_phi_deta3";

    data_stat->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_stat);
    if (delta_phi_deta3_stat == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3",delta_phi_deta3_model);
    if (delta_phi_deta3_model == 0) { cout << "merged_model_unc_delta_phi_deta3 not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta3",delta_phi_deta3_jer);
    if (delta_phi_deta3_jer == 0) { cout << "merged_jer_unc_delta_phi_deta3 not found!" << endl; return; }
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D(delta_phi_deta3_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta3, delta_phi_deta3_stat, delta_phi_deta3_model, delta_phi_deta3_jer, uncorr_unc, 3, output_path_plots, label_out + "delta_phi_deta3", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta3 norm distribution
    if (detail) { cout<<"Delta phi Deta3 Norm"<<endl; }

    TH1D *delta_phi_deta3_norm_stat = 0;
    TH1D *delta_phi_deta3_norm_model = 0;
    TH1D *delta_phi_deta3_norm_jer = 0;
    TString delta_phi_deta3_norm_name_out = label_out + "delta_phi_deta3_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta3_norm",delta_phi_deta3_norm_stat);
    if (delta_phi_deta3_norm_stat == 0) { cout << "ak5PF_delta_phi_deta3_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_norm",delta_phi_deta3_norm_model);
    if (delta_phi_deta3_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta3_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta3_norm",delta_phi_deta3_norm_jer);
    if (delta_phi_deta3_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta3_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_norm;
    delta_phi_deta3_norm =  new TH1D(delta_phi_deta3_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta3_norm, delta_phi_deta3_norm_stat, delta_phi_deta3_norm_model, delta_phi_deta3_norm_jer, uncorr_unc_norm, 3, output_path_plots, label_out + "delta_phi_deta3_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi Deta4"<<endl; }

    TH1D *delta_phi_deta4_stat = 0;
    TH1D *delta_phi_deta4_model = 0;
    TH1D *delta_phi_deta4_jer = 0;
    TString delta_phi_deta4_name_out = label_out + "delta_phi_deta4";

    data_stat->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_stat);
    if (delta_phi_deta4_stat == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4",delta_phi_deta4_model);
    if (delta_phi_deta4_model == 0) { cout << "merged_model_unc_delta_phi_deta4 not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta4",delta_phi_deta4_jer);
    if (delta_phi_deta4_jer == 0) { cout << "merged_jer_unc_delta_phi_deta4 not found!" << endl; return; }
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D(delta_phi_deta4_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta4, delta_phi_deta4_stat, delta_phi_deta4_model, delta_phi_deta4_jer, uncorr_unc, 4, output_path_plots, label_out + "delta_phi_deta4", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta4 norm distribution
    if (detail) { cout<<"Delta phi Deta4 Norm"<<endl; }

    TH1D *delta_phi_deta4_norm_stat = 0;
    TH1D *delta_phi_deta4_norm_model = 0;
    TH1D *delta_phi_deta4_norm_jer = 0;
    TString delta_phi_deta4_norm_name_out = label_out + "delta_phi_deta4_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta4_norm",delta_phi_deta4_norm_stat);
    if (delta_phi_deta4_norm_stat == 0) { cout << "ak5PF_delta_phi_deta4_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_norm",delta_phi_deta4_norm_model);
    if (delta_phi_deta4_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta4_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta4_norm",delta_phi_deta4_norm_jer);
    if (delta_phi_deta4_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta4_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_norm;
    delta_phi_deta4_norm =  new TH1D(delta_phi_deta4_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta4_norm, delta_phi_deta4_norm_stat, delta_phi_deta4_norm_model, delta_phi_deta4_norm_jer, uncorr_unc_norm, 4, output_path_plots, label_out + "delta_phi_deta4_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi Gap"<<endl; }

    TH1D *delta_phi_gap_stat = 0;
    TH1D *delta_phi_gap_model = 0;
    TH1D *delta_phi_gap_jer = 0;
    TString delta_phi_gap_name_out = label_out + "delta_phi_gap";

    data_stat->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_stat);
    if (delta_phi_gap_stat == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_gap",delta_phi_gap_model);
    if (delta_phi_gap_model == 0) { cout << "merged_model_unc_delta_phi_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_gap",delta_phi_gap_jer);
    if (delta_phi_gap_jer == 0) { cout << "merged_jer_unc_delta_phi_gap not found!" << endl; return; }
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D(delta_phi_gap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_gap, delta_phi_gap_stat, delta_phi_gap_model, delta_phi_gap_jer, uncorr_unc, 5, output_path_plots, label_out + "delta_phi_gap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi gap norm distribution
    if (detail) { cout<<"Delta phi Gap Norm"<<endl; }

    TH1D *delta_phi_gap_norm_stat = 0;
    TH1D *delta_phi_gap_norm_model = 0;
    TH1D *delta_phi_gap_norm_jer = 0;
    TString delta_phi_gap_norm_name_out = label_out + "delta_phi_gap_norm";

    data_stat->GetObject("ak5PF_delta_phi_gap_norm",delta_phi_gap_norm_stat);
    if (delta_phi_gap_norm_stat == 0) { cout << "ak5PF_delta_phi_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_gap_norm",delta_phi_gap_norm_model);
    if (delta_phi_gap_norm_model == 0) { cout << "merged_model_unc_delta_phi_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_gap_norm",delta_phi_gap_norm_jer);
    if (delta_phi_gap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_gap_norm;
    delta_phi_gap_norm =  new TH1D(delta_phi_gap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_gap_norm, delta_phi_gap_norm_stat, delta_phi_gap_norm_model, delta_phi_gap_norm_jer, uncorr_unc_norm, 5, output_path_plots, label_out + "delta_phi_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi Deta1 Gap"<<endl; }

    TH1D *delta_phi_deta1_gap_stat = 0;
    TH1D *delta_phi_deta1_gap_model = 0;
    TH1D *delta_phi_deta1_gap_jer = 0;
    TString delta_phi_deta1_gap_name_out = label_out + "delta_phi_deta1_gap";

    data_stat->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_stat);
    if (delta_phi_deta1_gap_stat == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_gap",delta_phi_deta1_gap_model);
    if (delta_phi_deta1_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta1_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta1_gap",delta_phi_deta1_gap_jer);
    if (delta_phi_deta1_gap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta1_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D(delta_phi_deta1_gap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta1_gap, delta_phi_deta1_gap_stat, delta_phi_deta1_gap_model, delta_phi_deta1_gap_jer, uncorr_unc, 6, output_path_plots, label_out + "delta_phi_deta1_gap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta1 gap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Gap Norm"<<endl; }

    TH1D *delta_phi_deta1_gap_norm_stat = 0;
    TH1D *delta_phi_deta1_gap_norm_model = 0;
    TH1D *delta_phi_deta1_gap_norm_jer = 0;
    TString delta_phi_deta1_gap_norm_name_out = label_out + "delta_phi_deta1_gap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_stat);
    if (delta_phi_deta1_gap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta1_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_model);
    if (delta_phi_deta1_gap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta1_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_jer);
    if (delta_phi_deta1_gap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta1_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap_norm;
    delta_phi_deta1_gap_norm =  new TH1D(delta_phi_deta1_gap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_stat, delta_phi_deta1_gap_norm_model, delta_phi_deta1_gap_norm_jer, uncorr_unc_norm, 6, output_path_plots, label_out + "delta_phi_deta1_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi Deta2 Gap"<<endl; }

    TH1D *delta_phi_deta2_gap_stat = 0;
    TH1D *delta_phi_deta2_gap_model = 0;
    TH1D *delta_phi_deta2_gap_jer = 0;
    TString delta_phi_deta2_gap_name_out = label_out + "delta_phi_deta2_gap";

    data_stat->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_stat);
    if (delta_phi_deta2_gap_stat == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_gap",delta_phi_deta2_gap_model);
    if (delta_phi_deta2_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta2_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta2_gap",delta_phi_deta2_gap_jer);
    if (delta_phi_deta2_gap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta2_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D(delta_phi_deta2_gap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta2_gap, delta_phi_deta2_gap_stat, delta_phi_deta2_gap_model, delta_phi_deta2_gap_jer, uncorr_unc, 7, output_path_plots, label_out + "delta_phi_deta2_gap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta2 gap norm stribution
    if (detail) { cout<<"Delta phi Deta2 Gap Norm"<<endl; }

    TH1D *delta_phi_deta2_gap_norm_stat = 0;
    TH1D *delta_phi_deta2_gap_norm_model = 0;
    TH1D *delta_phi_deta2_gap_norm_jer = 0;
    TString delta_phi_deta2_gap_norm_name_out = label_out + "delta_phi_deta2_gap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_stat);
    if (delta_phi_deta2_gap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta2_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_model);
    if (delta_phi_deta2_gap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta2_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_jer);
    if (delta_phi_deta2_gap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta2_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap_norm;
    delta_phi_deta2_gap_norm =  new TH1D(delta_phi_deta2_gap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_stat, delta_phi_deta2_gap_norm_model, delta_phi_deta2_gap_norm_jer, uncorr_unc_norm, 7, output_path_plots, label_out + "delta_phi_deta2_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi Deta3 Gap"<<endl; }

    TH1D *delta_phi_deta3_gap_stat = 0;
    TH1D *delta_phi_deta3_gap_model = 0;
    TH1D *delta_phi_deta3_gap_jer = 0;
    TString delta_phi_deta3_gap_name_out = label_out + "delta_phi_deta3_gap";

    data_stat->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_stat);
    if (delta_phi_deta3_gap_stat == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_gap",delta_phi_deta3_gap_model);
    if (delta_phi_deta3_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta3_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta3_gap",delta_phi_deta3_gap_jer);
    if (delta_phi_deta3_gap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta3_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D(delta_phi_deta3_gap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta3_gap, delta_phi_deta3_gap_stat, delta_phi_deta3_gap_model, delta_phi_deta3_gap_jer, uncorr_unc, 8, output_path_plots, label_out + "delta_phi_deta3_gap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta3 gap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Gap Norm"<<endl; }

    TH1D *delta_phi_deta3_gap_norm_stat = 0;
    TH1D *delta_phi_deta3_gap_norm_model = 0;
    TH1D *delta_phi_deta3_gap_norm_jer = 0;
    TString delta_phi_deta3_gap_norm_name_out = label_out + "delta_phi_deta3_gap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_stat);
    if (delta_phi_deta3_gap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta3_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_model);
    if (delta_phi_deta3_gap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta3_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_jer);
    if (delta_phi_deta3_gap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta3_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap_norm;
    delta_phi_deta3_gap_norm =  new TH1D(delta_phi_deta3_gap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_stat, delta_phi_deta3_gap_norm_model, delta_phi_deta3_gap_norm_jer, uncorr_unc_norm, 8, output_path_plots, label_out + "delta_phi_deta3_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi Deta4 Gap"<<endl; }

    TH1D *delta_phi_deta4_gap_stat = 0;
    TH1D *delta_phi_deta4_gap_model = 0;
    TH1D *delta_phi_deta4_gap_jer = 0;
    TString delta_phi_deta4_gap_name_out = label_out + "delta_phi_deta4_gap";

    data_stat->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_stat);
    if (delta_phi_deta4_gap_stat == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_gap",delta_phi_deta4_gap_model);
    if (delta_phi_deta4_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta4_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta4_gap",delta_phi_deta4_gap_jer);
    if (delta_phi_deta4_gap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta4_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D(delta_phi_deta4_gap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta4_gap, delta_phi_deta4_gap_stat, delta_phi_deta4_gap_model, delta_phi_deta4_gap_jer, uncorr_unc, 9, output_path_plots, label_out + "delta_phi_deta4_gap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta4 gap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Gap Norm"<<endl; }

    TH1D *delta_phi_deta4_gap_norm_stat = 0;
    TH1D *delta_phi_deta4_gap_norm_model = 0;
    TH1D *delta_phi_deta4_gap_norm_jer = 0;
    TString delta_phi_deta4_gap_norm_name_out = label_out + "delta_phi_deta4_gap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_stat);
    if (delta_phi_deta4_gap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta4_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_model);
    if (delta_phi_deta4_gap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta4_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_jer);
    if (delta_phi_deta4_gap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta4_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap_norm;
    delta_phi_deta4_gap_norm =  new TH1D(delta_phi_deta4_gap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_stat, delta_phi_deta4_gap_norm_model, delta_phi_deta4_gap_norm_jer, uncorr_unc_norm, 9, output_path_plots, label_out + "delta_phi_deta4_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi Nogap"<<endl; }

    TH1D *delta_phi_nogap_stat = 0;
    TH1D *delta_phi_nogap_model = 0;
    TH1D *delta_phi_nogap_jer = 0;
    TString delta_phi_nogap_name_out = label_out + "delta_phi_nogap";

    data_stat->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_stat);
    if (delta_phi_nogap_stat == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_nogap",delta_phi_nogap_model);
    if (delta_phi_nogap_model == 0) { cout << "merged_model_unc_delta_phi_nogap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_nogap",delta_phi_nogap_jer);
    if (delta_phi_nogap_jer == 0) { cout << "merged_jer_unc_delta_phi_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D(delta_phi_nogap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_nogap, delta_phi_nogap_stat, delta_phi_nogap_model, delta_phi_nogap_jer, uncorr_unc, 10, output_path_plots, label_out + "delta_phi_nogap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi nogap norm distribution
    if (detail) { cout<<"Delta phi Nogap Norm"<<endl; }

    TH1D *delta_phi_nogap_norm_stat = 0;
    TH1D *delta_phi_nogap_norm_model = 0;
    TH1D *delta_phi_nogap_norm_jer = 0;
    TString delta_phi_nogap_norm_name_out = label_out + "delta_phi_nogap_norm";

    data_stat->GetObject("ak5PF_delta_phi_nogap_norm",delta_phi_nogap_norm_stat);
    if (delta_phi_nogap_norm_stat == 0) { cout << "ak5PF_delta_phi_nogap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_nogap_norm",delta_phi_nogap_norm_model);
    if (delta_phi_nogap_norm_model == 0) { cout << "merged_model_unc_delta_phi_nogap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_nogap_norm",delta_phi_nogap_norm_jer);
    if (delta_phi_nogap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_nogap_norm;
    delta_phi_nogap_norm =  new TH1D(delta_phi_nogap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_nogap_norm, delta_phi_nogap_norm_stat, delta_phi_nogap_norm_model, delta_phi_nogap_norm_jer, uncorr_unc_norm, 10, output_path_plots, label_out + "delta_phi_nogap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_stat = 0;
    TH1D *delta_phi_deta1_nogap_model = 0;
    TH1D *delta_phi_deta1_nogap_jer = 0;
    TString delta_phi_deta1_nogap_name_out = label_out + "delta_phi_deta1_nogap";

    data_stat->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_stat);
    if (delta_phi_deta1_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_nogap",delta_phi_deta1_nogap_model);
    if (delta_phi_deta1_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta1_nogap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta1_nogap",delta_phi_deta1_nogap_jer);
    if (delta_phi_deta1_nogap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta1_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D(delta_phi_deta1_nogap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta1_nogap, delta_phi_deta1_nogap_stat, delta_phi_deta1_nogap_model, delta_phi_deta1_nogap_jer, uncorr_unc, 11, output_path_plots, label_out + "delta_phi_deta1_nogap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta1 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta1_nogap_norm_stat = 0;
    TH1D *delta_phi_deta1_nogap_norm_model = 0;
    TH1D *delta_phi_deta1_nogap_norm_jer = 0;
    TString delta_phi_deta1_nogap_norm_name_out = label_out + "delta_phi_deta1_nogap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_stat);
    if (delta_phi_deta1_nogap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_model);
    if (delta_phi_deta1_nogap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_jer);
    if (delta_phi_deta1_nogap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap_norm;
    delta_phi_deta1_nogap_norm =  new TH1D(delta_phi_deta1_nogap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_stat, delta_phi_deta1_nogap_norm_model, delta_phi_deta1_nogap_norm_jer, uncorr_unc_norm, 11, output_path_plots, label_out + "delta_phi_deta1_nogap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_stat = 0;
    TH1D *delta_phi_deta2_nogap_model = 0;
    TH1D *delta_phi_deta2_nogap_jer = 0;
    TString delta_phi_deta2_nogap_name_out = label_out + "delta_phi_deta2_nogap";

    data_stat->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_stat);
    if (delta_phi_deta2_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_nogap",delta_phi_deta2_nogap_model);
    if (delta_phi_deta2_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta2_nogap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta2_nogap",delta_phi_deta2_nogap_jer);
    if (delta_phi_deta2_nogap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta2_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D(delta_phi_deta2_nogap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta2_nogap, delta_phi_deta2_nogap_stat, delta_phi_deta2_nogap_model, delta_phi_deta2_nogap_jer, uncorr_unc, 12, output_path_plots, label_out + "delta_phi_deta2_nogap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta2 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta2_nogap_norm_stat = 0;
    TH1D *delta_phi_deta2_nogap_norm_model = 0;
    TH1D *delta_phi_deta2_nogap_norm_jer = 0;
    TString delta_phi_deta2_nogap_norm_name_out = label_out + "delta_phi_deta2_nogap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_stat);
    if (delta_phi_deta2_nogap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_model);
    if (delta_phi_deta2_nogap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_jer);
    if (delta_phi_deta2_nogap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap_norm;
    delta_phi_deta2_nogap_norm =  new TH1D(delta_phi_deta2_nogap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_stat, delta_phi_deta2_nogap_norm_model, delta_phi_deta2_nogap_norm_jer, uncorr_unc_norm, 12, output_path_plots, label_out + "delta_phi_deta2_nogap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_stat = 0;
    TH1D *delta_phi_deta3_nogap_model = 0;
    TH1D *delta_phi_deta3_nogap_jer = 0;
    TString delta_phi_deta3_nogap_name_out = label_out + "delta_phi_deta3_nogap";

    data_stat->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_stat);
    if (delta_phi_deta3_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_nogap",delta_phi_deta3_nogap_model);
    if (delta_phi_deta3_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta3_nogap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta3_nogap",delta_phi_deta3_nogap_jer);
    if (delta_phi_deta3_nogap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta3_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D(delta_phi_deta3_nogap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta3_nogap, delta_phi_deta3_nogap_stat, delta_phi_deta3_nogap_model, delta_phi_deta3_nogap_jer, uncorr_unc, 13, output_path_plots, label_out + "delta_phi_deta3_nogap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta3 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta3_nogap_norm_stat = 0;
    TH1D *delta_phi_deta3_nogap_norm_model = 0;
    TH1D *delta_phi_deta3_nogap_norm_jer = 0;
    TString delta_phi_deta3_nogap_norm_name_out = label_out + "delta_phi_deta3_nogap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_stat);
    if (delta_phi_deta3_nogap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_model);
    if (delta_phi_deta3_nogap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_jer);
    if (delta_phi_deta3_nogap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap_norm;
    delta_phi_deta3_nogap_norm =  new TH1D(delta_phi_deta3_nogap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_stat, delta_phi_deta3_nogap_norm_model, delta_phi_deta3_nogap_norm_jer, uncorr_unc_norm, 13, output_path_plots, label_out + "delta_phi_deta3_nogap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_stat = 0;
    TH1D *delta_phi_deta4_nogap_model = 0;
    TH1D *delta_phi_deta4_nogap_jer = 0;
    TString delta_phi_deta4_nogap_name_out = label_out + "delta_phi_deta4_nogap";

    data_stat->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_stat);
    if (delta_phi_deta4_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_nogap",delta_phi_deta4_nogap_model);
    if (delta_phi_deta4_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta4_nogap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta4_nogap",delta_phi_deta4_nogap_jer);
    if (delta_phi_deta4_nogap_jer == 0) { cout << "merged_jer_unc_delta_phi_deta4_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D(delta_phi_deta4_nogap_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta4_nogap, delta_phi_deta4_nogap_stat, delta_phi_deta4_nogap_model, delta_phi_deta4_nogap_jer, uncorr_unc, 14, output_path_plots, label_out + "delta_phi_deta4_nogap", "top_right", detail);


//compute uncorrelated uncertainty for delta phi deta4 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta4_nogap_norm_stat = 0;
    TH1D *delta_phi_deta4_nogap_norm_model = 0;
    TH1D *delta_phi_deta4_nogap_norm_jer = 0;
    TString delta_phi_deta4_nogap_norm_name_out = label_out + "delta_phi_deta4_nogap_norm";

    data_stat->GetObject("ak5PF_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_stat);
    if (delta_phi_deta4_nogap_norm_stat == 0) { cout << "ak5PF_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_model);
    if (delta_phi_deta4_nogap_norm_model == 0) { cout << "merged_model_unc_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_jer);
    if (delta_phi_deta4_nogap_norm_jer == 0) { cout << "merged_jer_unc_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap_norm;
    delta_phi_deta4_nogap_norm =  new TH1D(delta_phi_deta4_nogap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#phi [rad];Relative Uncorrelated Uncertainty", dphi_nbins, dphi_bins);

    calc_uncorr_uncertainty(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_stat, delta_phi_deta4_nogap_norm_model, delta_phi_deta4_nogap_norm_jer, uncorr_unc_norm, 14, output_path_plots, label_out + "delta_phi_deta4_nogap_norm", "top_right", detail);


//compute uncorrelated uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pT Inside Gap"<<endl; }

    TH1D *leading_pt_inside_gap_stat = 0;
    TH1D *leading_pt_inside_gap_model = 0;
    TH1D *leading_pt_inside_gap_jer = 0;
    TString leading_pt_inside_gap_name_out = label_out + "leading_pt_inside_gap";

    data_stat->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_stat);
    if (leading_pt_inside_gap_stat == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_inside_gap",leading_pt_inside_gap_model);
    if (leading_pt_inside_gap_model == 0) { cout << "merged_model_unc_leading_pt_inside_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_leading_pt_inside_gap",leading_pt_inside_gap_jer);
    if (leading_pt_inside_gap_jer == 0) { cout << "merged_jer_unc_leading_pt_inside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D(leading_pt_inside_gap_name_out,"Relative Uncorrelated Uncertainty;p_{T}^{inside} [GeV];Relative Uncorrelated Uncertainty", in_nbins, in_bins);

    calc_uncorr_uncertainty(leading_pt_inside_gap, leading_pt_inside_gap_stat, leading_pt_inside_gap_model, leading_pt_inside_gap_jer, uncorr_unc, 15, output_path_plots, label_out + "leading_pt_inside_gap", "top_right", detail);


//compute uncorrelated uncertainty for leading pt inside gap norm distribution
    if (detail) { cout<<"Leading pT Inside Gap Norm"<<endl; }

    TH1D *leading_pt_inside_gap_norm_stat = 0;
    TH1D *leading_pt_inside_gap_norm_model = 0;
    TH1D *leading_pt_inside_gap_norm_jer = 0;
    TString leading_pt_inside_gap_norm_name_out = label_out + "leading_pt_inside_gap_norm";

    data_stat->GetObject("ak5PF_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_stat);
    if (leading_pt_inside_gap_norm_stat == 0) { cout << "ak5PF_leading_pt_inside_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_model);
    if (leading_pt_inside_gap_norm_model == 0) { cout << "merged_model_unc_leading_pt_inside_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_jer);
    if (leading_pt_inside_gap_norm_jer == 0) { cout << "merged_jer_unc_leading_pt_inside_gap_norm not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap_norm;
    leading_pt_inside_gap_norm =  new TH1D(leading_pt_inside_gap_norm_name_out,"Relative Uncorrelated Uncertainty;p_{T}^{inside} [GeV];Relative Uncorrelated Uncertainty", in_nbins, in_bins);

    calc_uncorr_uncertainty(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_stat, leading_pt_inside_gap_norm_model, leading_pt_inside_gap_norm_jer, uncorr_unc_norm, 15, output_path_plots, label_out + "leading_pt_inside_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for leading eta star inside gap distribution
    if (detail) { cout<<"Leading Eta* Inside Gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_stat = 0;
    TH1D *leading_eta_star_inside_gap_model = 0;
    TH1D *leading_eta_star_inside_gap_jer = 0;
    TString leading_eta_star_inside_gap_name_out = label_out + "leading_eta_star_inside_gap";

    data_stat->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_stat);
    if (leading_eta_star_inside_gap_stat == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_eta_star_inside_gap",leading_eta_star_inside_gap_model);
    if (leading_eta_star_inside_gap_model == 0) { cout << "merged_model_unc_leading_eta_star_inside_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_leading_eta_star_inside_gap",leading_eta_star_inside_gap_jer);
    if (leading_eta_star_inside_gap_jer == 0) { cout << "merged_jer_unc_leading_eta_star_inside_gap not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D(leading_eta_star_inside_gap_name_out,"Relative Uncorrelated Uncertainty;#eta*;Relative Uncorrelated Uncertainty", etastar_nbins, etastar_bins);

    calc_uncorr_uncertainty(leading_eta_star_inside_gap, leading_eta_star_inside_gap_stat, leading_eta_star_inside_gap_model, leading_eta_star_inside_gap_jer, uncorr_unc, 16, output_path_plots, label_out + "leading_eta_star_inside_gap", "top_right", detail);


//compute uncorrelated uncertainty for leading eta star inside gap norm distribution
    if (detail) { cout<<"Leading Eta* Inside Gap Norm"<<endl; }

    TH1D *leading_eta_star_inside_gap_norm_stat = 0;
    TH1D *leading_eta_star_inside_gap_norm_model = 0;
    TH1D *leading_eta_star_inside_gap_norm_jer = 0;
    TString leading_eta_star_inside_gap_norm_name_out = label_out + "leading_eta_star_inside_gap_norm";

    data_stat->GetObject("ak5PF_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_stat);
    if (leading_eta_star_inside_gap_norm_stat == 0) { cout << "ak5PF_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_model);
    if (leading_eta_star_inside_gap_norm_model == 0) { cout << "merged_model_unc_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_jer);
    if (leading_eta_star_inside_gap_norm_jer == 0) { cout << "merged_jer_unc_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap_norm;
    leading_eta_star_inside_gap_norm =  new TH1D(leading_eta_star_inside_gap_norm_name_out,"Relative Uncorrelated Uncertainty;#eta*;Relative Uncorrelated Uncertainty", etastar_nbins, etastar_bins);

    calc_uncorr_uncertainty(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_stat, leading_eta_star_inside_gap_norm_model, leading_eta_star_inside_gap_norm_jer, uncorr_unc_norm, 16, output_path_plots, label_out + "leading_eta_star_inside_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Delta Eta Outside Gap"<<endl; }

    TH1D *delta_eta_outside_gap_stat = 0;
    TH1D *delta_eta_outside_gap_model = 0;
    TH1D *delta_eta_outside_gap_jer = 0;
    TString delta_eta_outside_gap_name_out = label_out + "delta_eta_outside_gap";

    data_stat->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_stat);
    if (delta_eta_outside_gap_stat == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_eta_outside_gap",delta_eta_outside_gap_model);
    if (delta_eta_outside_gap_model == 0) { cout << "merged_model_unc_delta_eta_outside_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_eta_outside_gap",delta_eta_outside_gap_jer);
    if (delta_eta_outside_gap_jer == 0) { cout << "merged_jer_unc_delta_eta_outside_gap not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D(delta_eta_outside_gap_name_out,"Relative Uncorrelated Uncertainty;#Delta#eta^{out};Relative Uncorrelated Uncertainty", deta_out_nbins, deta_out_bins);

    calc_uncorr_uncertainty(delta_eta_outside_gap, delta_eta_outside_gap_stat, delta_eta_outside_gap_model, delta_eta_outside_gap_jer, uncorr_unc, 17, output_path_plots, label_out + "delta_eta_outside_gap", "top_right", detail);


//compute uncorrelated uncertainty for delta eta outside gap norm distribution
    if (detail) { cout<<"Delta Eta Outside Gap Norm"<<endl; }

    TH1D *delta_eta_outside_gap_norm_stat = 0;
    TH1D *delta_eta_outside_gap_norm_model = 0;
    TH1D *delta_eta_outside_gap_norm_jer = 0;
    TString delta_eta_outside_gap_norm_name_out = label_out + "delta_eta_outside_gap_norm";

    data_stat->GetObject("ak5PF_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_stat);
    if (delta_eta_outside_gap_norm_stat == 0) { cout << "ak5PF_delta_eta_outside_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_model);
    if (delta_eta_outside_gap_norm_model == 0) { cout << "merged_model_unc_delta_eta_outside_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_jer);
    if (delta_eta_outside_gap_norm_jer == 0) { cout << "merged_jer_unc_delta_eta_outside_gap_norm not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap_norm;
    delta_eta_outside_gap_norm =  new TH1D(delta_eta_outside_gap_norm_name_out,"Relative Uncorrelated Uncertainty;#Delta#eta^{out};Relative Uncorrelated Uncertainty", deta_out_nbins, deta_out_bins);

    calc_uncorr_uncertainty(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_stat, delta_eta_outside_gap_norm_model, delta_eta_outside_gap_norm_jer, uncorr_unc_norm, 17, output_path_plots, label_out + "delta_eta_outside_gap_norm", "top_right", detail);


//compute uncorrelated uncertainty for leading pt outside gap distribution
    if (detail) { cout<<"Leading pT Outside Gap"<<endl; }

    TH1D *leading_pt_outside_gap_stat = 0;
    TH1D *leading_pt_outside_gap_model = 0;
    TH1D *leading_pt_outside_gap_jer = 0;
    TString leading_pt_outside_gap_name_out = label_out + "leading_pt_outside_gap";

    data_stat->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_stat);
    if (leading_pt_outside_gap_stat == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_outside_gap",leading_pt_outside_gap_model);
    if (leading_pt_outside_gap_model == 0) { cout << "merged_model_unc_leading_pt_outside_gap not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_leading_pt_outside_gap",leading_pt_outside_gap_jer);
    if (leading_pt_outside_gap_jer == 0) { cout << "merged_jer_unc_leading_pt_outside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D(leading_pt_outside_gap_name_out,"Relative Uncorrelated Uncertainty;p_{T}^{outside} [GeV];Relative Uncorrelated Uncertainty", out_nbins, out_bins);

    calc_uncorr_uncertainty(leading_pt_outside_gap, leading_pt_outside_gap_stat, leading_pt_outside_gap_model, leading_pt_outside_gap_jer, uncorr_unc, 18, output_path_plots, label_out + "leading_pt_outside_gap", "top_right", detail);


//compute uncorrelated uncertainty for leading pt outside gap norm distribution
    if (detail) { cout<<"Leading pT Outside Gap Norm"<<endl; }

    TH1D *leading_pt_outside_gap_norm_stat = 0;
    TH1D *leading_pt_outside_gap_norm_model = 0;
    TH1D *leading_pt_outside_gap_norm_jer = 0;
    TString leading_pt_outside_gap_norm_name_out = label_out + "leading_pt_outside_gap_norm";

    data_stat->GetObject("ak5PF_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_stat);
    if (leading_pt_outside_gap_norm_stat == 0) { cout << "ak5PF_leading_pt_outside_gap_norm not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_model);
    if (leading_pt_outside_gap_norm_model == 0) { cout << "merged_model_unc_leading_pt_outside_gap_norm not found!" << endl; return; }
    data_jer->GetObject("merged_jer_unc_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_jer);
    if (leading_pt_outside_gap_norm_jer == 0) { cout << "merged_jer_unc_leading_pt_outside_gap_norm not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap_norm;
    leading_pt_outside_gap_norm =  new TH1D(leading_pt_outside_gap_norm_name_out,"Relative Uncorrelated Uncertainty;p_{T}^{outside} [GeV];Relative Uncorrelated Uncertainty", out_nbins, out_bins);

    calc_uncorr_uncertainty(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_stat, leading_pt_outside_gap_norm_model, leading_pt_outside_gap_norm_jer, uncorr_unc_norm, 18, output_path_plots, label_out + "leading_pt_outside_gap_norm", "top_right", detail);


//output the error variation
    if (detail) { cout<<"Display Uncorrelated Uncertainties..."<<endl; }
    if (disp_uncertainty) { show_uncorr_uncertainties(uncorr_unc); }
    if (detail) { cout<<"Display Uncorrelated Uncertainties for Normalized Distributions..."<<endl; }
    if (disp_uncertainty) { show_uncorr_uncertainties(uncorr_unc_norm); }


//Opening the output root file
    if (detail) { cout<<"Creating " << uncorr_uncertainty << "..."<<endl; }
    TFile *data_output = TFile::Open( uncorr_uncertainty.c_str() , "RECREATE");


//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file..."<<endl; }
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
    leading_eta_star_inside_gap->Write();
    delta_eta_outside_gap->Write();
    leading_pt_outside_gap->Write();

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
    leading_eta_star_inside_gap_norm->Write();
    delta_eta_outside_gap_norm->Write();
    leading_pt_outside_gap_norm->Write();
    if (detail) { cout<<"Writing was sucessfull!"<<endl; }


//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data_stat->Close();
    data_model->Close();
    data_output->Close();


//deleting all histograms to avoid memory leak - causes a memory leak!
/*    delete(delta_phi);
    delete(delta_phi_deta1);
    delete(delta_phi_deta2);
    delete(delta_phi_deta3);
    delete(delta_phi_deta4);
    delete(delta_phi_gap);
    delete(delta_phi_deta1_gap);
    delete(delta_phi_deta2_gap);
    delete(delta_phi_deta3_gap);
    delete(delta_phi_deta4_gap);
    delete(delta_phi_nogap);
    delete(delta_phi_deta1_nogap);
    delete(delta_phi_deta2_nogap);
    delete(delta_phi_deta3_nogap);
    delete(delta_phi_deta4_nogap);
    delete(leading_pt_inside_gap);
    delete(leading_pt_outside_gap);
    delete(delta_eta_outside_gap);
    delete(leading_eta_star_inside_gap); */

    if (detail) { cout<<"Done!"<<endl; }


}
