// Pedro Cipriano, Nov 2013
// DESY, CMS
// Last Update: 06 Nov 2013
//
// compute_corr_uncertainty()

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


void calc_corr_uncertainty(TH1 *corr, TH1 *jes, double *corr_unc, int index = 0, string path = "../output/" , string fileout = "test", string legend_position = "top_left", bool detail = false)
{
//initialize variables
double max = 0.0;
double min = 0.0;
double tot = 0.0;
double ave = 0.0;

//set base values for the uncertainties
double val_lumi = 0.04;
double val_trigger = 0.01;
double val_jes = 0.0;
double val_corr = 0.0;

//setting the mode
//if (mode == "up") { val_trigger = 0.01;}
//if (mode == "down") { val_trigger = 0.00;}

//clonning the histograms
TH1D *lumi =  (TH1D*) jes->Clone();
TH1D *trigger = (TH1D*) jes->Clone();

//loop over the bins
for (int i = 1; i <= corr->GetNbinsX();i++)
	{
	val_jes = jes->GetBinContent(i);
	val_corr = sqrt(val_jes * val_jes +  val_lumi * val_lumi + val_trigger * val_trigger);
	corr->SetBinContent(i,val_corr);
	lumi->SetBinContent(i,val_lumi);
	trigger->SetBinContent(i,val_trigger);
        if (val_corr > max) {max = val_corr;}
        tot = tot + val_corr;
        if (val_corr < min || i == 1) { min = val_corr;}
        if (detail) { cout<< "Bin " << i << " Jes = " << val_jes << ", Lumi = " << val_lumi << ", Trigger = " << val_trigger << " and Correlated total = " << val_corr << endl; }
	}

//calculates the average model uncertainty
    ave = tot/corr->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    corr_unc[index*3+0] = ave*100;
    corr_unc[index*3+1] = min*100;
    corr_unc[index*3+2] = max*100;

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
    double max_plot = 1.0;

//plooting
    jes->SetMaximum(max_plot);
    corr->SetMaximum(max_plot);
    corr->SetMinimum(min_plot);
    format_histogram(corr, 1, 1);
    corr->Draw("hist");
    format_histogram(jes, 2, 2);
    jes->Draw("hist same");
    format_histogram(lumi, 4, 3);
    lumi->Draw("hist same");
    format_histogram(trigger, 6, 4);
    trigger->Draw("hist same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    int curves = 4;
    set_legend_position(legend_position, curves, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(corr,"Correlated Uncertainty","l");
    leg01->AddEntry(jes,"Jet Energy Scale Uncertainty","l");
    leg01->AddEntry(lumi,"Luminosity Uncertainty","l");
    leg01->AddEntry(trigger,"Trigger Uncertainty","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

//print the plot
    print_plots(c01, path, fileout);

}


void calc_corr_uncertainty_norm(TH1 *corr, TH1 *jes, double *corr_unc, int index = 0, string path = "../output/" , string fileout = "test", string legend_position = "top_left", bool detail = false)
{
//initialize variables
double max = 0.0;
double min = 0.0;
double tot = 0.0;
double ave = 0.0;

//set base values for the uncertainties
double val_trigger = 0.01;
double val_jes = 0.0;
double val_corr = 0.0;

//setting the mode
//if (mode == "up") { val_trigger = 0.01;}
//if (mode == "down") { val_trigger = 0.01;}

//clonning the histograms
TH1D *trigger = (TH1D*) jes->Clone();

//loop over the bins
for (int i = 1; i <= corr->GetNbinsX();i++)
	{
	val_jes = jes->GetBinContent(i);
	val_corr = sqrt(val_jes * val_jes + val_trigger * val_trigger);
	corr->SetBinContent(i,val_corr);
	trigger->SetBinContent(i,val_trigger);
        if (val_corr > max) {max = val_corr;}
        tot = tot + val_corr;
        if (val_corr < min || i == 1) { min = val_corr;}
        if (detail) { cout<< "Bin " << i << " Jes = " << val_jes << ", Trigger = " << val_trigger << " and Correlated total = " << val_corr << endl; }
	}

//calculates the average model uncertainty
    ave = tot/corr->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    corr_unc[index*3+0] = ave*100;
    corr_unc[index*3+1] = min*100;
    corr_unc[index*3+2] = max*100;

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
    double max_plot = 1.0;

//plooting
    jes->SetMaximum(max_plot);
    corr->SetMaximum(max_plot);
    corr->SetMinimum(min_plot);
    format_histogram(corr, 1, 1);
    corr->Draw("hist");
    format_histogram(jes, 2, 2);
    jes->Draw("hist same");
    format_histogram(trigger, 6, 4);
    trigger->Draw("hist same");


//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    int curves = 3;
    set_legend_position(legend_position, curves, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(corr,"Correlated Uncertainty","l");
    leg01->AddEntry(jes,"Jet Energy Scale Uncertainty","l");
    leg01->AddEntry(trigger,"Trigger Uncertainty","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

//print the plot
    print_plots(c01, path, fileout);

}

void show_corr_uncertainties(double *corr_unc)
{
//shows the computed jes uncertainties
    cout<<" "<<endl;
    cout<<"Total Uncertainty"<<endl;
    cout<<"Observable                  Average  Minimum  Maximum"<<endl;
    cout<<"Delta phi                   "<<corr_unc[0]<<"  "<<corr_unc[1]<<" "<<corr_unc[2]<<endl;
    cout<<"Delta phi deta1             "<<corr_unc[3]<<"  "<<corr_unc[4]<<"  "<<corr_unc[5]<<endl;
    cout<<"Delta phi deta2             "<<corr_unc[6]<<"   "<<corr_unc[7]<<" "<<corr_unc[8]<<endl;
    cout<<"Delta phi deta3             "<<corr_unc[9]<<"  "<<corr_unc[10]<<" "<<corr_unc[11]<<endl;
    cout<<"Delta phi deta4             "<<corr_unc[12]<<"   "<<corr_unc[13]<<" "<<corr_unc[14]<<endl;
    cout<<"Delta phi gap               "<<corr_unc[15]<<"  "<<corr_unc[16]<<" "<<corr_unc[17]<<endl;
    cout<<"Delta phi deta1 gap         "<<corr_unc[18]<<"  "<<corr_unc[19]<<"  "<<corr_unc[20]<<endl;
    cout<<"Delta phi deta2 gap         "<<corr_unc[21]<<"  "<<corr_unc[22]<<" "<<corr_unc[23]<<endl;
    cout<<"Delta phi deta3 gap         "<<corr_unc[24]<<"  "<<corr_unc[25]<<" "<<corr_unc[26]<<endl;
    cout<<"Delta phi deta4 gap         "<<corr_unc[27]<<"  "<<corr_unc[28]<<" "<<corr_unc[29]<<endl;
    cout<<"Delta phi nogap             "<<corr_unc[30]<<"  "<<corr_unc[31]<<" "<<corr_unc[32]<<endl;
    cout<<"Delta phi deta1 nogap       "<<corr_unc[33]<<"  "<<corr_unc[34]<<"  "<<corr_unc[35]<<endl;
    cout<<"Delta phi deta2 nogap       "<<corr_unc[36]<<"   "<<corr_unc[37]<<"  "<<corr_unc[38]<<endl;
    cout<<"Delta phi deta3 nogap       "<<corr_unc[39]<<"  "<<corr_unc[40]<<"  "<<corr_unc[41]<<endl;
    cout<<"Delta phi deta4 nogap       "<<corr_unc[42]<<"  "<<corr_unc[43]<<"  "<<corr_unc[44]<<endl;
    cout<<"Leading pT inside gap       "<<corr_unc[45]<<"  "<<corr_unc[46]<<" "<<corr_unc[47]<<endl;
    cout<<"Leading eta star inside gap "<<corr_unc[48]<<"  "<<corr_unc[49]<<"  "<<corr_unc[50]<<endl;
    cout<<"Delta eta outside gap       "<<corr_unc[51]<<"  "<<corr_unc[52]<<" "<<corr_unc[53]<<endl;
    cout<<"Leading pT outside gap      "<<corr_unc[54]<<"  "<<corr_unc[55]<<" "<<corr_unc[56]<<endl;
}


void corr_uncertainty(string path_jes, string label_in, string corr_uncertainty, string label_out, string output_path_plots, bool detail = false, bool disp_uncertainty = true, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Total Uncertainty Configuration"<<endl; }
//    if (detail) { cout << "Mode:                 " << mode << endl; }
    if (detail) { cout << "Input path for JES:   " << path_jes << endl; }
    if (detail) { cout << "Label In:             " << label_in << endl; }
    if (detail) { cout << "Label Out:            " << label_out << endl; }
    if (detail) { cout << "Output path:          " << corr_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:    " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:         " << detail << endl; }
    if (detail) { cout << "Display Results:      " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:            " << test << endl; }

//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_jes = new TFile( path_jes.c_str() );

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
    double corr_unc[19*3];
    double corr_unc_norm[19*3];

    for (int i=0; i<= 19*3-1;i++)
    {
    corr_unc[i] = 0.0;
    corr_unc_norm[i] = 0.0;
    }


//compute total uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_jes = 0;
    TString delta_phi_name_in = label_in + "delta_phi";
    TString delta_phi_name_out = label_out + "delta_phi";

    data_jes->GetObject(delta_phi_name_in,delta_phi_jes);
    if (delta_phi_jes == 0) { cout << delta_phi_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi;
    delta_phi =  new TH1D(delta_phi_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi, delta_phi_jes, corr_unc, 0, output_path_plots, label_out + "delta_phi", "top_right", detail);

//compute total uncertainty for delta phi norm distribution
    if (detail) { cout<<"Delta phi Norm"<<endl; }

    TH1D *delta_phi_norm_jes = 0;
    TString delta_phi_norm_name_in = label_in + "delta_phi_norm";
    TString delta_phi_norm_name_out = label_out + "delta_phi_norm";

    data_jes->GetObject(delta_phi_norm_name_in,delta_phi_norm_jes);
    if (delta_phi_norm_jes == 0) { cout << delta_phi_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_norm;
    delta_phi_norm =  new TH1D(delta_phi_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_norm, delta_phi_norm_jes, corr_unc_norm, 0, output_path_plots, label_out + "delta_phi_norm", "top_right", detail);


//compute total uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi Deta1"<<endl; }

    TH1D *delta_phi_deta1_jes = 0;
    TString delta_phi_deta1_name_in = label_in + "delta_phi_deta1";
    TString delta_phi_deta1_name_out = label_out + "delta_phi_deta1";

    data_jes->GetObject(delta_phi_deta1_name_in,delta_phi_deta1_jes);
    if (delta_phi_deta1_jes == 0) { cout << delta_phi_deta1_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D(delta_phi_deta1_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta1, delta_phi_deta1_jes, corr_unc, 1, output_path_plots, label_out + "delta_phi_deta1", "top_right", detail);


//compute total uncertainty for delta phi deta1 norm distribution
    if (detail) { cout<<"Delta phi Deta1 Norm"<<endl; }

    TH1D *delta_phi_deta1_norm_jes = 0;
    TString delta_phi_deta1_norm_name_in = label_in + "delta_phi_deta1_norm";
    TString delta_phi_deta1_norm_name_out = label_out + "delta_phi_deta1_norm";

    data_jes->GetObject(delta_phi_deta1_norm_name_in,delta_phi_deta1_norm_jes);
    if (delta_phi_deta1_norm_jes == 0) { cout << delta_phi_deta1_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_norm;
    delta_phi_deta1_norm =  new TH1D(delta_phi_deta1_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta1_norm, delta_phi_deta1_norm_jes, corr_unc_norm, 1, output_path_plots, label_out + "delta_phi_deta1_norm", "top_right", detail);



//compute total uncertainty for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi Deta2"<<endl; }

    TH1D *delta_phi_deta2_jes = 0;
    TString delta_phi_deta2_name_in = label_in + "delta_phi_deta2";
    TString delta_phi_deta2_name_out = label_out + "delta_phi_deta2";

    data_jes->GetObject(delta_phi_deta2_name_in,delta_phi_deta2_jes);
    if (delta_phi_deta2_jes == 0) { cout << delta_phi_deta2_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D(delta_phi_deta2_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta2, delta_phi_deta2_jes, corr_unc, 2, output_path_plots, label_out + "delta_phi_deta2", "top_right", detail);


//compute total uncertainty for delta phi deta2 norm distribution
    if (detail) { cout<<"Delta phi Deta2 Norm"<<endl; }

    TH1D *delta_phi_deta2_norm_jes = 0;
    TString delta_phi_deta2_norm_name_in = label_in + "delta_phi_deta2_norm";
    TString delta_phi_deta2_norm_name_out = label_out + "delta_phi_deta2_norm";

    data_jes->GetObject(delta_phi_deta2_norm_name_in,delta_phi_deta2_norm_jes);
    if (delta_phi_deta2_norm_jes == 0) { cout << delta_phi_deta2_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_norm;
    delta_phi_deta2_norm =  new TH1D(delta_phi_deta2_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta2_norm, delta_phi_deta2_norm_jes, corr_unc_norm, 2, output_path_plots, label_out + "delta_phi_deta2_norm", "top_right", detail);


//compute total uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi Deta3"<<endl; }

    TH1D *delta_phi_deta3_jes = 0;
    TString delta_phi_deta3_name_in = label_in + "delta_phi_deta3";
    TString delta_phi_deta3_name_out = label_out + "delta_phi_deta3";

    data_jes->GetObject(delta_phi_deta3_name_in,delta_phi_deta3_jes);
    if (delta_phi_deta3_jes == 0) { cout << delta_phi_deta3_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D(delta_phi_deta3_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta3, delta_phi_deta3_jes, corr_unc, 3, output_path_plots, label_out + "delta_phi_deta3", "top_right", detail);


//compute total uncertainty for delta phi deta3 norm distribution
    if (detail) { cout<<"Delta phi Deta3 Norm"<<endl; }

    TH1D *delta_phi_deta3_norm_jes = 0;
    TString delta_phi_deta3_norm_name_in = label_in + "delta_phi_deta3_norm";
    TString delta_phi_deta3_norm_name_out = label_out + "delta_phi_deta3_norm";

    data_jes->GetObject(delta_phi_deta3_norm_name_in,delta_phi_deta3_norm_jes);
    if (delta_phi_deta3_norm_jes == 0) { cout << delta_phi_deta3_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_norm;
    delta_phi_deta3_norm =  new TH1D(delta_phi_deta3_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta3_norm, delta_phi_deta3_norm_jes, corr_unc_norm, 3, output_path_plots, label_out + "delta_phi_deta3_norm", "top_right", detail);


//compute total uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi Deta4"<<endl; }

    TH1D *delta_phi_deta4_jes = 0;
    TString delta_phi_deta4_name_in = label_in + "delta_phi_deta4";
    TString delta_phi_deta4_name_out = label_out + "delta_phi_deta4";

    data_jes->GetObject(delta_phi_deta4_name_in,delta_phi_deta4_jes);
    if (delta_phi_deta4_jes == 0) { cout << delta_phi_deta4_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D(delta_phi_deta4_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta4, delta_phi_deta4_jes, corr_unc, 4, output_path_plots, label_out + "delta_phi_deta4", "top_right", detail);


//compute total uncertainty for delta phi deta4 norm distribution
    if (detail) { cout<<"Delta phi Deta4 Norm"<<endl; }

    TH1D *delta_phi_deta4_norm_jes = 0;
    TString delta_phi_deta4_norm_name_in = label_in + "delta_phi_deta4_norm";
    TString delta_phi_deta4_norm_name_out = label_out + "delta_phi_deta4_norm";

    data_jes->GetObject(delta_phi_deta4_norm_name_in,delta_phi_deta4_norm_jes);
    if (delta_phi_deta4_norm_jes == 0) { cout << delta_phi_deta4_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_norm;
    delta_phi_deta4_norm =  new TH1D(delta_phi_deta4_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta4_norm, delta_phi_deta4_norm_jes, corr_unc_norm, 4, output_path_plots, label_out + "delta_phi_deta4_norm", "top_right", detail);


//compute total uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi Gap"<<endl; }

    TH1D *delta_phi_gap_jes = 0;
    TString delta_phi_gap_name_in = label_in + "delta_phi_gap";
    TString delta_phi_gap_name_out = label_out + "delta_phi_gap";

    data_jes->GetObject(delta_phi_gap_name_in,delta_phi_gap_jes);
    if (delta_phi_gap_jes == 0) { cout << delta_phi_gap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D(delta_phi_gap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_gap, delta_phi_gap_jes, corr_unc, 5, output_path_plots, label_out + "delta_phi_gap", "top_right", detail);


//compute total uncertainty for delta phi gap norm distribution
    if (detail) { cout<<"Delta phi Gap Norm"<<endl; }

    TH1D *delta_phi_gap_norm_jes = 0;
    TString delta_phi_gap_norm_name_in = label_in + "delta_phi_gap_norm";
    TString delta_phi_gap_norm_name_out = label_out + "delta_phi_gap_norm";

    data_jes->GetObject(delta_phi_gap_norm_name_in,delta_phi_gap_norm_jes);
    if (delta_phi_gap_norm_jes == 0) { cout << delta_phi_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_gap_norm;
    delta_phi_gap_norm =  new TH1D(delta_phi_gap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_gap_norm, delta_phi_gap_norm_jes, corr_unc_norm, 5, output_path_plots, label_out + "delta_phi_gap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi Deta1 Gap"<<endl; }

    TH1D *delta_phi_deta1_gap_jes = 0;
    TString delta_phi_deta1_gap_name_in = label_in + "delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_name_out = label_out + "delta_phi_deta1_gap";

    data_jes->GetObject(delta_phi_deta1_gap_name_in,delta_phi_deta1_gap_jes);
    if (delta_phi_deta1_gap_jes == 0) { cout << delta_phi_deta1_gap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D(delta_phi_deta1_gap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta1_gap, delta_phi_deta1_gap_jes, corr_unc, 6, output_path_plots, label_out + "delta_phi_deta1_gap", "top_right", detail);


//compute total uncertainty for delta phi deta1 gap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Gap Norm"<<endl; }

    TH1D *delta_phi_deta1_gap_norm_jes = 0;
    TString delta_phi_deta1_gap_norm_name_in = label_in + "delta_phi_deta1_gap_norm";
    TString delta_phi_deta1_gap_norm_name_out = label_out + "delta_phi_deta1_gap_norm";

    data_jes->GetObject(delta_phi_deta1_gap_norm_name_in,delta_phi_deta1_gap_norm_jes);
    if (delta_phi_deta1_gap_norm_jes == 0) { cout << delta_phi_deta1_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap_norm;
    delta_phi_deta1_gap_norm =  new TH1D(delta_phi_deta1_gap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_jes, corr_unc_norm, 6, output_path_plots, label_out + "delta_phi_deta1_gap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi Deta2 Gap"<<endl; }

    TH1D *delta_phi_deta2_gap_jes = 0;
    TString delta_phi_deta2_gap_name_in = label_in + "delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_name_out = label_out + "delta_phi_deta2_gap";

    data_jes->GetObject(delta_phi_deta2_gap_name_in,delta_phi_deta2_gap_jes);
    if (delta_phi_deta2_gap_jes == 0) { cout << delta_phi_deta2_gap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D(delta_phi_deta2_gap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta2_gap, delta_phi_deta2_gap_jes, corr_unc, 7, output_path_plots, label_out + "delta_phi_deta2_gap", "top_right", detail);


//compute total uncertainty for delta phi deta2 gap norm distribution
    if (detail) { cout<<"Delta phi Deta2 Gap Norm"<<endl; }

    TH1D *delta_phi_deta2_gap_norm_jes = 0;
    TString delta_phi_deta2_gap_norm_name_in = label_in + "delta_phi_deta2_gap_norm";
    TString delta_phi_deta2_gap_norm_name_out = label_out + "delta_phi_deta2_gap_norm";

    data_jes->GetObject(delta_phi_deta2_gap_norm_name_in,delta_phi_deta2_gap_norm_jes);
    if (delta_phi_deta2_gap_norm_jes == 0) { cout << delta_phi_deta2_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap_norm;
    delta_phi_deta2_gap_norm =  new TH1D(delta_phi_deta2_gap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_jes, corr_unc_norm, 7, output_path_plots, label_out + "delta_phi_deta2_gap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi Deta3 Gap"<<endl; }

    TH1D *delta_phi_deta3_gap_jes = 0;
    TString delta_phi_deta3_gap_name_in = label_in + "delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_name_out = label_out + "delta_phi_deta3_gap";

    data_jes->GetObject(delta_phi_deta3_gap_name_in,delta_phi_deta3_gap_jes);
    if (delta_phi_deta3_gap_jes == 0) { cout << delta_phi_deta3_gap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D(delta_phi_deta3_gap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta3_gap, delta_phi_deta3_gap_jes, corr_unc, 8, output_path_plots, label_out + "delta_phi_deta3_gap", "top_right", detail);


//compute total uncertainty for delta phi deta3 gap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Gap Norm"<<endl; }

    TH1D *delta_phi_deta3_gap_norm_jes = 0;
    TString delta_phi_deta3_gap_norm_name_in = label_in + "delta_phi_deta3_gap_norm";
    TString delta_phi_deta3_gap_norm_name_out = label_out + "delta_phi_deta3_gap_norm";

    data_jes->GetObject(delta_phi_deta3_gap_norm_name_in,delta_phi_deta3_gap_norm_jes);
    if (delta_phi_deta3_gap_norm_jes == 0) { cout << delta_phi_deta3_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap_norm;
    delta_phi_deta3_gap_norm =  new TH1D(delta_phi_deta3_gap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_jes, corr_unc_norm, 8, output_path_plots, label_out + "delta_phi_deta3_gap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi Deta4 Gap"<<endl; }

    TH1D *delta_phi_deta4_gap_jes = 0;
    TString delta_phi_deta4_gap_name_in = label_in + "delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_name_out = label_out + "delta_phi_deta4_gap";

    data_jes->GetObject(delta_phi_deta4_gap_name_in,delta_phi_deta4_gap_jes);
    if (delta_phi_deta4_gap_jes == 0) { cout << delta_phi_deta4_gap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D(delta_phi_deta4_gap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta4_gap, delta_phi_deta4_gap_jes, corr_unc, 9, output_path_plots, label_out + "delta_phi_deta4_gap", "top_right", detail);


//compute total uncertainty for delta phi deta4 gap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Gap Norm"<<endl; }

    TH1D *delta_phi_deta4_gap_norm_jes = 0;
    TString delta_phi_deta4_gap_norm_name_in = label_in + "delta_phi_deta4_gap_norm";
    TString delta_phi_deta4_gap_norm_name_out = label_out + "delta_phi_deta4_gap_norm";

    data_jes->GetObject(delta_phi_deta4_gap_norm_name_in,delta_phi_deta4_gap_norm_jes);
    if (delta_phi_deta4_gap_norm_jes == 0) { cout << delta_phi_deta4_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap_norm;
    delta_phi_deta4_gap_norm =  new TH1D(delta_phi_deta4_gap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_jes, corr_unc_norm, 9, output_path_plots, label_out + "delta_phi_deta4_gap_norm", "top_right", detail);


//compute total uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi Nogap"<<endl; }

    TH1D *delta_phi_nogap_jes = 0;
    TString delta_phi_nogap_name_in = label_in + "delta_phi_nogap";
    TString delta_phi_nogap_name_out = label_out + "delta_phi_nogap";

    data_jes->GetObject(delta_phi_nogap_name_in,delta_phi_nogap_jes);
    if (delta_phi_nogap_jes == 0) { cout << delta_phi_nogap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D(delta_phi_nogap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_nogap, delta_phi_nogap_jes, corr_unc, 10, output_path_plots, label_out + "delta_phi_nogap", "top_right", detail);


//compute total uncertainty for delta phi nogap norm distribution
    if (detail) { cout<<"Delta phi Nogap Norm"<<endl; }

    TH1D *delta_phi_nogap_norm_jes = 0;
    TString delta_phi_nogap_norm_name_in = label_in + "delta_phi_nogap_norm";
    TString delta_phi_nogap_norm_name_out = label_out + "delta_phi_nogap_norm";

    data_jes->GetObject(delta_phi_nogap_norm_name_in,delta_phi_nogap_norm_jes);
    if (delta_phi_nogap_norm_jes == 0) { cout << delta_phi_nogap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_nogap_norm;
    delta_phi_nogap_norm =  new TH1D(delta_phi_nogap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_nogap_norm, delta_phi_nogap_norm_jes, corr_unc_norm, 10, output_path_plots, label_out + "delta_phi_nogap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_jes = 0;
    TString delta_phi_deta1_nogap_name_in = label_in + "delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_name_out = label_out + "delta_phi_deta1_nogap";

    data_jes->GetObject(delta_phi_deta1_nogap_name_in,delta_phi_deta1_nogap_jes);
    if (delta_phi_deta1_nogap_jes == 0) { cout << delta_phi_deta1_nogap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D(delta_phi_deta1_nogap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jes, corr_unc, 11, output_path_plots, label_out + "delta_phi_deta1_nogap", "top_right", detail);


//compute total uncertainty for delta phi deta1 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta1_nogap_norm_jes = 0;
    TString delta_phi_deta1_nogap_norm_name_in = label_in + "delta_phi_deta1_nogap_norm";
    TString delta_phi_deta1_nogap_norm_name_out = label_out + "delta_phi_deta1_nogap_norm";

    data_jes->GetObject(delta_phi_deta1_nogap_norm_name_in,delta_phi_deta1_nogap_norm_jes);
    if (delta_phi_deta1_nogap_norm_jes == 0) { cout << delta_phi_deta1_nogap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap_norm;
    delta_phi_deta1_nogap_norm =  new TH1D(delta_phi_deta1_nogap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_jes, corr_unc_norm, 11, output_path_plots, label_out + "delta_phi_deta1_nogap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_jes = 0;
    TString delta_phi_deta2_nogap_name_in = label_in + "delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_name_out = label_out + "delta_phi_deta2_nogap";

    data_jes->GetObject(delta_phi_deta2_nogap_name_in,delta_phi_deta2_nogap_jes);
    if (delta_phi_deta2_nogap_jes == 0) { cout << delta_phi_deta2_nogap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D(delta_phi_deta2_nogap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jes, corr_unc, 12, output_path_plots, label_out + "delta_phi_deta2_nogap", "top_right", detail);


//compute total uncertainty for delta phi deta2 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta2_nogap_norm_jes = 0;
    TString delta_phi_deta2_nogap_norm_name_in = label_in + "delta_phi_deta2_nogap_norm";
    TString delta_phi_deta2_nogap_norm_name_out = label_out + "delta_phi_deta2_nogap_norm";

    data_jes->GetObject(delta_phi_deta2_nogap_norm_name_in,delta_phi_deta2_nogap_norm_jes);
    if (delta_phi_deta2_nogap_norm_jes == 0) { cout << delta_phi_deta2_nogap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap_norm;
    delta_phi_deta2_nogap_norm =  new TH1D(delta_phi_deta2_nogap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_jes, corr_unc_norm, 12, output_path_plots, label_out + "delta_phi_deta2_nogap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_jes = 0;
    TString delta_phi_deta3_nogap_name_in = label_in + "delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_name_out = label_out + "delta_phi_deta3_nogap";

    data_jes->GetObject(delta_phi_deta3_nogap_name_in,delta_phi_deta3_nogap_jes);
    if (delta_phi_deta3_nogap_jes == 0) { cout << delta_phi_deta3_nogap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D(delta_phi_deta3_nogap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jes, corr_unc, 13, output_path_plots, label_out + "delta_phi_deta3_nogap", "top_right", detail);


//compute total uncertainty for delta phi deta3 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta3_nogap_norm_jes = 0;
    TString delta_phi_deta3_nogap_norm_name_in = label_in + "delta_phi_deta3_nogap_norm";
    TString delta_phi_deta3_nogap_norm_name_out = label_out + "delta_phi_deta3_nogap_norm";

    data_jes->GetObject(delta_phi_deta3_nogap_norm_name_in,delta_phi_deta3_nogap_norm_jes);
    if (delta_phi_deta3_nogap_norm_jes == 0) { cout << delta_phi_deta3_nogap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap_norm;
    delta_phi_deta3_nogap_norm =  new TH1D(delta_phi_deta3_nogap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_jes, corr_unc_norm, 13, output_path_plots, label_out + "delta_phi_deta3_nogap_norm", "top_right", detail);


//compute total uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_jes = 0;
    TString delta_phi_deta4_nogap_name_in = label_in + "delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_name_out = label_out + "delta_phi_deta4_nogap";

    data_jes->GetObject(delta_phi_deta4_nogap_name_in,delta_phi_deta4_nogap_jes);
    if (delta_phi_deta4_nogap_jes == 0) { cout << delta_phi_deta4_nogap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D(delta_phi_deta4_nogap_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jes, corr_unc, 14, output_path_plots, label_out + "delta_phi_deta4_nogap", "top_right", detail);


//compute total uncertainty for delta phi deta4 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta4_nogap_norm_jes = 0;
    TString delta_phi_deta4_nogap_norm_name_in = label_in + "delta_phi_deta4_nogap_norm";
    TString delta_phi_deta4_nogap_norm_name_out = label_out + "delta_phi_deta4_nogap_norm";

    data_jes->GetObject(delta_phi_deta4_nogap_norm_name_in,delta_phi_deta4_nogap_norm_jes);
    if (delta_phi_deta4_nogap_norm_jes == 0) { cout << delta_phi_deta4_nogap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap_norm;
    delta_phi_deta4_nogap_norm =  new TH1D(delta_phi_deta4_nogap_norm_name_out,"Relative Correlated Uncertainty;#Delta#phi [rad];Relative Correlated Uncertainty", dphi_nbins, dphi_bins);

    calc_corr_uncertainty_norm(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_jes, corr_unc_norm, 14, output_path_plots, label_out + "delta_phi_deta4_nogap_norm", "top_right", detail);



//compute total uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pT Inside Gap"<<endl; }

    TH1D *leading_pt_inside_gap_jes = 0;
    TString leading_pt_inside_gap_name_in = label_in + "leading_pt_inside_gap";
    TString leading_pt_inside_gap_name_out = label_out + "leading_pt_inside_gap";

    data_jes->GetObject(leading_pt_inside_gap_name_in,leading_pt_inside_gap_jes);
    if (leading_pt_inside_gap_jes == 0) { cout << leading_pt_inside_gap_name_in << " not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D(leading_pt_inside_gap_name_out,"Relative Correlated Uncertainty;p_{T}^{inside} [GeV];Relative Correlated Uncertainty", in_nbins, in_bins);

    calc_corr_uncertainty(leading_pt_inside_gap, leading_pt_inside_gap_jes, corr_unc, 15, output_path_plots, label_out + "leading_pt_inside_gap", "top_right", detail);


//compute total uncertainty for leading pt inside gap norm distribution
    if (detail) { cout<<"Leading pT Inside Gap Norm"<<endl; }

    TH1D *leading_pt_inside_gap_norm_jes = 0;
    TString leading_pt_inside_gap_norm_name_in = label_in + "leading_pt_inside_gap_norm";
    TString leading_pt_inside_gap_norm_name_out = label_out + "leading_pt_inside_gap_norm";

    data_jes->GetObject(leading_pt_inside_gap_norm_name_in,leading_pt_inside_gap_norm_jes);
    if (leading_pt_inside_gap_norm_jes == 0) { cout << leading_pt_inside_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap_norm;
    leading_pt_inside_gap_norm =  new TH1D(leading_pt_inside_gap_norm_name_out,"Relative Correlated Uncertainty;p_{T}^{inside} [GeV];Relative Correlated Uncertainty", in_nbins, in_bins);

    calc_corr_uncertainty(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_jes, corr_unc_norm, 15, output_path_plots, label_out + "leading_pt_inside_gap_norm", "top_right", detail);


//compute total uncertainty for leading eta star inside gap distribution
    if (detail) { cout<<"Leading Eta* Inside Gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_jes = 0;
    TString leading_eta_star_inside_gap_name_in = label_in + "leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_name_out = label_out + "leading_eta_star_inside_gap";

    data_jes->GetObject(leading_eta_star_inside_gap_name_in,leading_eta_star_inside_gap_jes);
    if (leading_eta_star_inside_gap_jes == 0) { cout << leading_eta_star_inside_gap_name_in << " not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D(leading_eta_star_inside_gap_name_out,"Relative Correlated Uncertainty;#eta*;Relative Correlated Uncertainty", etastar_nbins, etastar_bins);

    calc_corr_uncertainty(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jes, corr_unc, 16, output_path_plots, label_out + "leading_eta_star_inside_gap", "top_right", detail);


//compute total uncertainty for leading eta star inside gap norm distribution
    if (detail) { cout<<"Leading Eta* Inside Gap Norm"<<endl; }

    TH1D *leading_eta_star_inside_gap_norm_jes = 0;
    TString leading_eta_star_inside_gap_norm_name_in = label_in + "leading_eta_star_inside_gap_norm";
    TString leading_eta_star_inside_gap_norm_name_out = label_out + "leading_eta_star_inside_gap_norm";

    data_jes->GetObject(leading_eta_star_inside_gap_norm_name_in,leading_eta_star_inside_gap_norm_jes);
    if (leading_eta_star_inside_gap_norm_jes == 0) { cout << leading_eta_star_inside_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap_norm;
    leading_eta_star_inside_gap_norm =  new TH1D(leading_eta_star_inside_gap_norm_name_out,"Relative Correlated Uncertainty;#eta*;Relative Correlated Uncertainty", etastar_nbins, etastar_bins);

    calc_corr_uncertainty_norm(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_jes, corr_unc_norm, 16, output_path_plots, label_out + "leading_eta_star_inside_gap_norm", "top_right", detail);


//compute total uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Delta Eta Outside Gap"<<endl; }

    TH1D *delta_eta_outside_gap_jes = 0;
    TString delta_eta_outside_gap_name_in = label_in + "delta_eta_outside_gap";
    TString delta_eta_outside_gap_name_out = label_out + "delta_eta_outside_gap";

    data_jes->GetObject(delta_eta_outside_gap_name_in,delta_eta_outside_gap_jes);
    if (delta_eta_outside_gap_jes == 0) { cout << delta_eta_outside_gap_name_in << " not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D(delta_eta_outside_gap_name_out,"Relative Correlated Uncertainty;#Delta#eta^{out};Relative Correlated Uncertainty", deta_out_nbins, deta_out_bins);

    calc_corr_uncertainty(delta_eta_outside_gap, delta_eta_outside_gap_jes, corr_unc, 17, output_path_plots, label_out + "delta_eta_outside_gap", "top_right", detail);


//compute total uncertainty for delta eta outside gap norm distribution
    if (detail) { cout<<"Delta Eta Outside Gap Norm"<<endl; }

    TH1D *delta_eta_outside_gap_norm_jes = 0;
    TString delta_eta_outside_gap_norm_name_in = label_in + "delta_eta_outside_gap_norm";
    TString delta_eta_outside_gap_norm_name_out = label_out + "delta_eta_outside_gap_norm";

    data_jes->GetObject(delta_eta_outside_gap_norm_name_in,delta_eta_outside_gap_norm_jes);
    if (delta_eta_outside_gap_norm_jes == 0) { cout << delta_eta_outside_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap_norm;
    delta_eta_outside_gap_norm =  new TH1D(delta_eta_outside_gap_norm_name_out,"Relative Correlated Uncertainty;#Delta#eta^{out};Relative Correlated Uncertainty", deta_out_nbins, deta_out_bins);

    calc_corr_uncertainty_norm(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_jes, corr_unc_norm, 17, output_path_plots, label_out + "delta_eta_outside_gap_norm", "top_right", detail);


//compute total uncertainty for leading pt outside gap distribution
    if (detail) { cout<<"Leading pT Outside Gap"<<endl; }

    TH1D *leading_pt_outside_gap_jes = 0;
    TString leading_pt_outside_gap_name_in = label_in + "leading_pt_outside_gap";
    TString leading_pt_outside_gap_name_out = label_out + "leading_pt_outside_gap";

    data_jes->GetObject(leading_pt_outside_gap_name_in,leading_pt_outside_gap_jes);
    if (leading_pt_outside_gap_jes == 0) { cout << leading_pt_outside_gap_name_in << " not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D(leading_pt_outside_gap_name_out,"Relative Correlated Uncertainty;p_{T}^{outside} [GeV];Relative Correlated Uncertainty", out_nbins, out_bins);

    calc_corr_uncertainty(leading_pt_outside_gap, leading_pt_outside_gap_jes, corr_unc, 18, output_path_plots, label_out + "leading_pt_outside_gap", "top_right", detail);


//compute total uncertainty for leading pt outside gap norm distribution
    if (detail) { cout<<"Leading pT Outside Gap Norm"<<endl; }

    TH1D *leading_pt_outside_gap_norm_jes = 0;
    TString leading_pt_outside_gap_norm_name_in = label_in + "leading_pt_outside_gap_norm";
    TString leading_pt_outside_gap_norm_name_out = label_out + "leading_pt_outside_gap_norm";

    data_jes->GetObject(leading_pt_outside_gap_norm_name_in,leading_pt_outside_gap_norm_jes);
    if (leading_pt_outside_gap_norm_jes == 0) { cout << leading_pt_outside_gap_norm_name_in << " not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap_norm;
    leading_pt_outside_gap_norm =  new TH1D(leading_pt_outside_gap_norm_name_out,"Relative Correlated Uncertainty;p_{T}^{outside} [GeV];Relative Correlated Uncertainty", out_nbins, out_bins);

    calc_corr_uncertainty_norm(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_jes, corr_unc_norm, 18, output_path_plots, label_out + "leading_pt_outside_gap_norm", "top_right", detail);


//output the error variation
    if (detail) { cout<<"Display Correlated Uncertainties..."<<endl; }
    if (disp_uncertainty) { show_corr_uncertainties(corr_unc); }
    if (detail) { cout<<"Display Correlated Uncertainties for Normalized Distributions..."<<endl; }
    if (disp_uncertainty) { show_corr_uncertainties(corr_unc_norm); }


//Opening the output root file
    if (detail) { cout<<"Creating " << corr_uncertainty << "..."<<endl; }
    TFile *data_output = TFile::Open( corr_uncertainty.c_str() , "RECREATE");


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
    data_jes->Close();
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
