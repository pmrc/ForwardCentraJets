// Pedro Cipriano, Mar 2011
// DESY, CMS
// Last Update: 22 Mar 2013
//
// total_uncertainty()

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

void calc_total_uncertainty(TH1D *total, TH1D *jes, TH1D *model, TH1D *data, double *total_unc, int index, bool detail)
{

double max = 0.0;
double min = 0.0;
double tot = 0.0;
double ave = 0.0;

double val_lumi = 0.04;
double val_jes = 0.0;
double val_model = 0.0;
double val_stat = 0.0;
double val_total = 0.0;

for (int i = 1; i <= data->GetNbinsX();i++)
	{
	val_jes = jes->GetBinContent(i);
	val_model = model->GetBinContent(i);
	val_stat = data->GetBinError(i)/data->GetBinContent(i);
	val_total = sqrt(val_jes * val_jes + val_model * val_model + val_stat * val_stat + val_lumi * val_lumi);
	total->SetBinContent(i,val_total);
        if (val_total > max) {max = val_total;}
        tot = tot + val_total;
        if (val_total < min || i == 1) { min = val_total;}
        if (detail) { cout<< "Bin " << i << " Jes = " << val_jes << ", Model = " << val_model << ", Stat = " << val_stat << " and Total = " << val_total << endl; }
	}

//calculates the average model uncertainty
    ave = tot/total->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    total_unc[index*3+0] = ave*100;
    total_unc[index*3+1] = min*100;
    total_unc[index*3+2] = max*100;

}


void show_total_uncertainties(double *total_unc)
{
//shows the computed jes uncertainties
    cout<<" "<<endl;
    cout<<"Total Uncertainty"<<endl;
    cout<<"Observable                  Average  Minimum  Maximum"<<endl;
    cout<<"Delta phi                   "<<total_unc[0]<<"  "<<total_unc[1]<<" "<<total_unc[2]<<endl;
    cout<<"Delta phi deta1             "<<total_unc[3]<<"  "<<total_unc[4]<<"  "<<total_unc[5]<<endl;
    cout<<"Delta phi deta2             "<<total_unc[6]<<"   "<<total_unc[7]<<" "<<total_unc[8]<<endl;
    cout<<"Delta phi deta3             "<<total_unc[9]<<"  "<<total_unc[10]<<" "<<total_unc[11]<<endl;
    cout<<"Delta phi deta4             "<<total_unc[12]<<"   "<<total_unc[13]<<" "<<total_unc[14]<<endl;
    cout<<"Delta phi gap               "<<total_unc[15]<<"  "<<total_unc[16]<<" "<<total_unc[17]<<endl;
    cout<<"Delta phi deta1 gap         "<<total_unc[18]<<"  "<<total_unc[19]<<"  "<<total_unc[20]<<endl;
    cout<<"Delta phi deta2 gap         "<<total_unc[21]<<"  "<<total_unc[22]<<" "<<total_unc[23]<<endl;
    cout<<"Delta phi deta3 gap         "<<total_unc[24]<<"  "<<total_unc[25]<<" "<<total_unc[26]<<endl;
    cout<<"Delta phi deta4 gap         "<<total_unc[27]<<"  "<<total_unc[28]<<" "<<total_unc[29]<<endl;
    cout<<"Delta phi nogap             "<<total_unc[30]<<"  "<<total_unc[31]<<" "<<total_unc[32]<<endl;
    cout<<"Delta phi deta1 nogap       "<<total_unc[33]<<"  "<<total_unc[34]<<"  "<<total_unc[35]<<endl;
    cout<<"Delta phi deta2 nogap       "<<total_unc[36]<<"   "<<total_unc[37]<<"  "<<total_unc[38]<<endl;
    cout<<"Delta phi deta3 nogap       "<<total_unc[39]<<"  "<<total_unc[40]<<"  "<<total_unc[41]<<endl;
    cout<<"Delta phi deta4 nogap       "<<total_unc[42]<<"  "<<total_unc[43]<<"  "<<total_unc[44]<<endl;
    cout<<"Leading pT inside gap       "<<total_unc[45]<<"  "<<total_unc[46]<<" "<<total_unc[47]<<endl;
    cout<<"Leading eta star inside gap "<<total_unc[48]<<"  "<<total_unc[49]<<"  "<<total_unc[50]<<endl;
    cout<<"Delta eta outside gap       "<<total_unc[51]<<"  "<<total_unc[52]<<" "<<total_unc[53]<<endl;
    cout<<"Leading pT outside gap      "<<total_unc[54]<<"  "<<total_unc[55]<<" "<<total_unc[56]<<endl;
}


void calc_total_uncertainty_v2(TH1D *total, TH1D *corr_up, TH1D *corr_down, TH1D *uncorr, double *total_unc_up, double *total_unc_down, int index, string legend_position, string path, string fileout, bool detail)
{

double max_up = 0.0, max_down = 0.0;
double min_up = 0.0, min_down = 0.0;
double tot_up = 0.0, tot_down = 0.0;
double ave_up = 0.0, ave_down = 0.0;

double val_corr_up, val_corr_down, val_uncorr, val_total_up, val_total_down;
double val_shift, val_unc;

//clone histograms
TH1D *corr_plot =  (TH1D*) total->Clone();
TH1D *uncorr_plot =  (TH1D*) total->Clone();

for (int i = 1; i <= total->GetNbinsX();i++)
	{
	val_corr_up = corr_up->GetBinContent(i);
	val_corr_down = corr_down->GetBinContent(i);
	val_uncorr = uncorr->GetBinContent(i);
	val_total_up = sqrt(val_corr_up * val_corr_up + val_uncorr * val_uncorr);
	val_total_down = sqrt(val_corr_down * val_corr_down + val_uncorr * val_uncorr);
        if (val_total_up > max_up) {max_up = val_total_up;}
        tot_up = tot_up + val_total_up;
        if (val_total_up < min_up || i == 1) { min_up = val_total_up;}
        if (val_total_down > max_down) {max_down = val_total_down;}
        tot_down = tot_down + val_total_down;
        if (val_total_down < min_down || i == 1) { min_down = val_total_down;}
        if (detail) { cout<< "Bin " << i << " Correlated Up = " << val_corr_up << ", Correlated Down = " << val_corr_down << ", Uncorrelated = " << val_uncorr << ", Total Up = " << val_total_up << " and Total Down = " << val_total_down <<endl; }
	uncorr_plot->SetBinContent(i,1);
	uncorr_plot->SetBinError(i,val_uncorr);
	val_unc = (val_corr_up + val_corr_down)/2.0;
	val_shift = (val_corr_up - val_corr_down)/2.0;
	corr_plot->SetBinContent(i,1+val_shift);
	corr_plot->SetBinError(i,val_unc);
	val_unc = (val_total_up + val_total_down)/2.0;
	val_shift = (val_total_up - val_total_down)/2.0;
	total->SetBinContent(i,1+val_shift);
	total->SetBinError(i,val_unc);
	}

//calculates the average model uncertainty
    ave_up = tot_up/total->GetNbinsX();
    ave_down = tot_down/total->GetNbinsX();

//displays the result
    if (detail)
	{
	cout<<"Result Up  : min ="<<min_up*100<<", max = "<<max_up*100<<", ave= "<<ave_up*100<<endl;
	cout<<"Result Down: min ="<<min_down*100<<", max = "<<max_down*100<<", ave= "<<ave_down*100<<endl;
	}

//saves the results in the results array
    total_unc_up[index*3+0] = ave_up*100;
    total_unc_up[index*3+1] = min_up*100;
    total_unc_up[index*3+2] = max_up*100;
    total_unc_down[index*3+0] = ave_down*100;
    total_unc_down[index*3+1] = min_down*100;
    total_unc_down[index*3+2] = max_down*100;


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


//plooting
    total->SetMaximum(2.7);
    total->SetMinimum(0.0);
    total->SetFillColor(6);
    total->Draw("e2");
    corr_plot->SetFillColor(5);
    corr_plot->Draw("e2 same");
    format_histogram(uncorr_plot, 1, 1);
    uncorr_plot->Draw("e1 same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 3, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(total,"Total Uncertainty","f");
    leg01->AddEntry(corr_plot,"Correlated Uncertainty","f");
    leg01->AddEntry(uncorr_plot,"Uncorrelated Uncertainty","lep");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);

}

void total_uncertainty_v2(string path_corr_up, string path_corr_down, string path_uncorr, string total_uncertainty, string output_path_plots, bool detail = false, bool disp_uncertainty = true, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Total Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for Correlated Uncertainties Up:   " << path_corr_up << endl; }
    if (detail) { cout << "Input path for Correlated Uncertainties Down: " << path_corr_down << endl; }
    if (detail) { cout << "Input path for Uncorrelated Uncertanties:     " << path_uncorr << endl; }
    if (detail) { cout << "Output path:                                  " << total_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:                            " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:                                 " << detail << endl; }
    if (detail) { cout << "Display Results:                              " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:                                    " << test << endl; }

//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_corr_up = new TFile( path_corr_up.c_str() );
    TFile *data_corr_down = new TFile( path_corr_down.c_str() );
    TFile *data_uncorr = new TFile( path_uncorr.c_str() );


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
    double total_unc_up[19*3], total_unc_down[19*3], total_unc_norm_up[19*3], total_unc_norm_down[19*3];

    for (int i=0; i<= 19*3-1;i++)
    	{
    	total_unc_up[i] = 0.0;
    	total_unc_down[i] = 0.0;
    	total_unc_norm_up[i] = 0.0;
    	total_unc_norm_down[i] = 0.0;
    	}


//compute total uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_corr_up = 0;
    TH1D *delta_phi_corr_down = 0;
    TH1D *delta_phi_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi",delta_phi_corr_up);
    if (delta_phi_corr_up == 0) { cout << "corr_unc_up_delta_phi not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi",delta_phi_corr_down);
    if (delta_phi_corr_down == 0) { cout << "corr_unc_down_delta_phi not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi",delta_phi_uncorr);
    if (delta_phi_uncorr == 0) { cout << "uncorr_unc_ddelta_phi not found!" << endl; return; }
    
    TH1D *delta_phi;
    delta_phi =  new TH1D("total_unc_delta_phi","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi, delta_phi_corr_up, delta_phi_corr_down, delta_phi_uncorr, total_unc_up, total_unc_down,  0, "top_right", output_path_plots, "tot_unc_delta_phi", detail);


//compute total uncertainty for delta phi norm distribution
    if (detail) { cout<<"Delta phi Norm"<<endl; }

    TH1D *delta_phi_norm_corr_up = 0;
    TH1D *delta_phi_norm_corr_down = 0;
    TH1D *delta_phi_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_norm",delta_phi_norm_corr_up);
    if (delta_phi_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_norm",delta_phi_norm_corr_down);
    if (delta_phi_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_norm",delta_phi_norm_uncorr);
    if (delta_phi_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_norm not found!" << endl; return; }
    
    TH1D *delta_phi_norm;
    delta_phi_norm =  new TH1D("total_unc_delta_phi_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_norm, delta_phi_norm_corr_up, delta_phi_norm_corr_down, delta_phi_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  0, "top_right", output_path_plots, "tot_unc_delta_phi_norm", detail);


//compute total uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi Deta1"<<endl; }

    TH1D *delta_phi_deta1_corr_up = 0;
    TH1D *delta_phi_deta1_corr_down = 0;
    TH1D *delta_phi_deta1_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta1",delta_phi_deta1_corr_up);
    if (delta_phi_deta1_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta1 not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta1",delta_phi_deta1_corr_down);
    if (delta_phi_deta1_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta1 not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta1",delta_phi_deta1_uncorr);
    if (delta_phi_deta1_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta1 not found!" << endl; return; }
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D("total_unc_delta_phi_deta1","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta1, delta_phi_deta1_corr_up, delta_phi_deta1_corr_down, delta_phi_deta1_uncorr, total_unc_up, total_unc_down,  1, "top_right", output_path_plots, "tot_unc_delta_phi_deta1", detail);


//compute total uncertainty for delta phi deta1 norm distribution
    if (detail) { cout<<"Delta phi Deta1 Norm"<<endl; }

    TH1D *delta_phi_deta1_norm_corr_up = 0;
    TH1D *delta_phi_deta1_norm_corr_down = 0;
    TH1D *delta_phi_deta1_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta1_norm",delta_phi_deta1_norm_corr_up);
    if (delta_phi_deta1_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta1_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta1_norm",delta_phi_deta1_norm_corr_down);
    if (delta_phi_deta1_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta1_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta1_norm",delta_phi_deta1_norm_uncorr);
    if (delta_phi_deta1_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta1_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_norm;
    delta_phi_deta1_norm =  new TH1D("total_unc_delta_phi_deta1_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta1_norm, delta_phi_deta1_norm_corr_up, delta_phi_deta1_norm_corr_down, delta_phi_deta1_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  1, "top_right", output_path_plots, "tot_unc_delta_phi_deta1_norm", detail);


//compute total uncertainty for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi Deta2"<<endl; }

    TH1D *delta_phi_deta2_corr_up = 0;
    TH1D *delta_phi_deta2_corr_down = 0;
    TH1D *delta_phi_deta2_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta2",delta_phi_deta2_corr_up);
    if (delta_phi_deta2_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta2 not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta2",delta_phi_deta2_corr_down);
    if (delta_phi_deta2_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta2 not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta2",delta_phi_deta2_uncorr);
    if (delta_phi_deta2_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta2 not found!" << endl; return; }
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D("total_unc_delta_phi_deta2","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta2, delta_phi_deta2_corr_up, delta_phi_deta2_corr_down, delta_phi_deta2_uncorr, total_unc_up, total_unc_down,  2, "top_right", output_path_plots, "tot_unc_delta_phi_deta2", detail);


//compute total uncertainty for delta phi deta2 norm distribution
    if (detail) { cout<<"Delta phi Deta2 Norm"<<endl; }

    TH1D *delta_phi_deta2_norm_corr_up = 0;
    TH1D *delta_phi_deta2_norm_corr_down = 0;
    TH1D *delta_phi_deta2_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta2_norm",delta_phi_deta2_norm_corr_up);
    if (delta_phi_deta2_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta2_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta2_norm",delta_phi_deta2_norm_corr_down);
    if (delta_phi_deta2_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta2_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta2_norm",delta_phi_deta2_norm_uncorr);
    if (delta_phi_deta2_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta2_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_norm;
    delta_phi_deta2_norm =  new TH1D("total_unc_delta_phi_deta2_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta2_norm, delta_phi_deta2_norm_corr_up, delta_phi_deta2_norm_corr_down, delta_phi_deta2_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  2, "top_right", output_path_plots, "tot_unc_delta_phi_deta2_norm", detail);


//compute total uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi Deta3"<<endl; }

    TH1D *delta_phi_deta3_corr_up = 0;
    TH1D *delta_phi_deta3_corr_down = 0;
    TH1D *delta_phi_deta3_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta3",delta_phi_deta3_corr_up);
    if (delta_phi_deta3_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta3 not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta3",delta_phi_deta3_corr_down);
    if (delta_phi_deta3_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta3 not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta3",delta_phi_deta3_uncorr);
    if (delta_phi_deta3_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta3 not found!" << endl; return; }
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D("total_unc_delta_phi_deta3","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta3, delta_phi_deta3_corr_up, delta_phi_deta3_corr_down, delta_phi_deta3_uncorr, total_unc_up, total_unc_down,  3, "top_right", output_path_plots, "tot_unc_delta_phi_deta3", detail);



//compute total uncertainty for delta phi deta3 norm distribution
    if (detail) { cout<<"Delta phi Deta3 Norm"<<endl; }

    TH1D *delta_phi_deta3_norm_corr_up = 0;
    TH1D *delta_phi_deta3_norm_corr_down = 0;
    TH1D *delta_phi_deta3_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta3_norm",delta_phi_deta3_norm_corr_up);
    if (delta_phi_deta3_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta3_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta3_norm",delta_phi_deta3_norm_corr_down);
    if (delta_phi_deta3_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta3_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta3_norm",delta_phi_deta3_norm_uncorr);
    if (delta_phi_deta3_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta3_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_norm;
    delta_phi_deta3_norm =  new TH1D("total_unc_delta_phi_deta3_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta3_norm, delta_phi_deta3_norm_corr_up, delta_phi_deta3_norm_corr_down, delta_phi_deta3_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  3, "top_right", output_path_plots, "tot_unc_delta_phi_deta3_norm", detail);


//compute total uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi Deta4"<<endl; }

    TH1D *delta_phi_deta4_corr_up = 0;
    TH1D *delta_phi_deta4_corr_down = 0;
    TH1D *delta_phi_deta4_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta4",delta_phi_deta4_corr_up);
    if (delta_phi_deta4_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta4 not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta4",delta_phi_deta4_corr_down);
    if (delta_phi_deta4_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta4 not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta4",delta_phi_deta4_uncorr);
    if (delta_phi_deta4_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta4 not found!" << endl; return; }
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D("total_unc_delta_phi_deta4","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta4, delta_phi_deta4_corr_up, delta_phi_deta4_corr_down, delta_phi_deta4_uncorr, total_unc_up, total_unc_down,  4, "top_right", output_path_plots, "tot_unc_delta_phi_deta4", detail);


//compute total uncertainty for delta phi deta4 norm distribution
    if (detail) { cout<<"Delta phi Deta4 Norm"<<endl; }

    TH1D *delta_phi_deta4_norm_corr_up = 0;
    TH1D *delta_phi_deta4_norm_corr_down = 0;
    TH1D *delta_phi_deta4_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta4_norm",delta_phi_deta4_norm_corr_up);
    if (delta_phi_deta4_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta4_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta4_norm",delta_phi_deta4_norm_corr_down);
    if (delta_phi_deta4_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta4_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta4_norm",delta_phi_deta4_norm_uncorr);
    if (delta_phi_deta4_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta4_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_norm;
    delta_phi_deta4_norm =  new TH1D("total_unc_delta_phi_deta4_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta4_norm, delta_phi_deta4_norm_corr_up, delta_phi_deta4_norm_corr_down, delta_phi_deta4_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  4, "top_right", output_path_plots, "tot_unc_delta_phi_deta4_norm", detail);


//compute total uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi Gap"<<endl; }

    TH1D *delta_phi_gap_corr_up = 0;
    TH1D *delta_phi_gap_corr_down = 0;
    TH1D *delta_phi_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_gap",delta_phi_gap_corr_up);
    if (delta_phi_gap_corr_up == 0) { cout << "corr_unc_up_delta_phi_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_gap",delta_phi_gap_corr_down);
    if (delta_phi_gap_corr_down == 0) { cout << "corr_unc_down_delta_phi_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_gap",delta_phi_gap_uncorr);
    if (delta_phi_gap_uncorr == 0) { cout << "uncorr_unc_delta_phi_gap not found!" << endl; return; }
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D("total_unc_delta_phi_gap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_gap, delta_phi_gap_corr_up, delta_phi_gap_corr_down, delta_phi_gap_uncorr, total_unc_up, total_unc_down,  5, "top_right", output_path_plots, "tot_unc_delta_phi_gap", detail);


//compute total uncertainty for delta phi gap norm distribution
    if (detail) { cout<<"Delta phi Gap Norm"<<endl; }

    TH1D *delta_phi_gap_norm_corr_up = 0;
    TH1D *delta_phi_gap_norm_corr_down = 0;
    TH1D *delta_phi_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_gap_norm",delta_phi_gap_norm_corr_up);
    if (delta_phi_gap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_gap_norm",delta_phi_gap_norm_corr_down);
    if (delta_phi_gap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_gap_norm",delta_phi_gap_norm_uncorr);
    if (delta_phi_gap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_gap_norm;
    delta_phi_gap_norm =  new TH1D("total_unc_delta_phi_gap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_gap_norm, delta_phi_gap_norm_corr_up, delta_phi_gap_norm_corr_down, delta_phi_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  5, "top_right", output_path_plots, "tot_unc_delta_phi_gap_norm", detail);


//compute total uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi Deta1 Gap"<<endl; }

    TH1D *delta_phi_deta1_gap_corr_up = 0;
    TH1D *delta_phi_deta1_gap_corr_down = 0;
    TH1D *delta_phi_deta1_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta1_gap",delta_phi_deta1_gap_corr_up);
    if (delta_phi_deta1_gap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta1_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta1_gap",delta_phi_deta1_gap_corr_down);
    if (delta_phi_deta1_gap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta1_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta1_gap",delta_phi_deta1_gap_uncorr);
    if (delta_phi_deta1_gap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta1_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D("total_unc_delta_phi_deta1_gap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta1_gap, delta_phi_deta1_gap_corr_up, delta_phi_deta1_gap_corr_down, delta_phi_deta1_gap_uncorr, total_unc_up, total_unc_down,  6, "top_right", output_path_plots, "tot_unc_delta_phi_deta1_gap", detail);


//compute total uncertainty for delta phi deta1 gap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Gap Norm"<<endl; }

    TH1D *delta_phi_deta1_gap_norm_corr_up = 0;
    TH1D *delta_phi_deta1_gap_norm_corr_down = 0;
    TH1D *delta_phi_deta1_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_corr_up);
    if (delta_phi_deta1_gap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta1_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_corr_down);
    if (delta_phi_deta1_gap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta1_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_uncorr);
    if (delta_phi_deta1_gap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta1_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap_norm;
    delta_phi_deta1_gap_norm =  new TH1D("total_unc_delta_phi_deta1_gap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_corr_up, delta_phi_deta1_gap_norm_corr_down, delta_phi_deta1_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  6, "top_right", output_path_plots, "tot_unc_delta_phi_deta1_gap_norm", detail);


//compute total uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi Deta2 Gap"<<endl; }

    TH1D *delta_phi_deta2_gap_corr_up = 0;
    TH1D *delta_phi_deta2_gap_corr_down = 0;
    TH1D *delta_phi_deta2_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta2_gap",delta_phi_deta2_gap_corr_up);
    if (delta_phi_deta2_gap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta2_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta2_gap",delta_phi_deta2_gap_corr_down);
    if (delta_phi_deta2_gap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta2_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta2_gap",delta_phi_deta2_gap_uncorr);
    if (delta_phi_deta2_gap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta2_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D("total_unc_delta_phi_deta2_gap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta2_gap, delta_phi_deta2_gap_corr_up, delta_phi_deta2_gap_corr_down, delta_phi_deta2_gap_uncorr, total_unc_up, total_unc_down,  7, "top_right", output_path_plots, "tot_unc_delta_phi_deta2_gap", detail);


//compute total uncertainty for delta phi deta2 gap norm distribution
    if (detail) { cout<<"Delta phi Deta2 Gap Norm"<<endl; }

    TH1D *delta_phi_deta2_gap_norm_corr_up = 0;
    TH1D *delta_phi_deta2_gap_norm_corr_down = 0;
    TH1D *delta_phi_deta2_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_corr_up);
    if (delta_phi_deta2_gap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta2_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_corr_down);
    if (delta_phi_deta2_gap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta2_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_uncorr);
    if (delta_phi_deta2_gap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta2_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap_norm;
    delta_phi_deta2_gap_norm =  new TH1D("total_unc_delta_phi_deta2_gap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_corr_up, delta_phi_deta2_gap_norm_corr_down, delta_phi_deta2_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  7, "top_right", output_path_plots, "tot_unc_delta_phi_deta2_gap_norm", detail);


//compute total uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi Deta3 Gap"<<endl; }

    TH1D *delta_phi_deta3_gap_corr_up = 0;
    TH1D *delta_phi_deta3_gap_corr_down = 0;
    TH1D *delta_phi_deta3_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta3_gap",delta_phi_deta3_gap_corr_up);
    if (delta_phi_deta3_gap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta3_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta3_gap",delta_phi_deta3_gap_corr_down);
    if (delta_phi_deta3_gap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta3_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta3_gap",delta_phi_deta3_gap_uncorr);
    if (delta_phi_deta3_gap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta3_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D("total_unc_delta_phi_deta3_gap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta3_gap, delta_phi_deta3_gap_corr_up, delta_phi_deta3_gap_corr_down, delta_phi_deta3_gap_uncorr, total_unc_up, total_unc_down,  8, "top_right", output_path_plots, "tot_unc_delta_phi_deta3_gap", detail);


//compute total uncertainty for delta phi deta3 gap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Gap Norm"<<endl; }

    TH1D *delta_phi_deta3_gap_norm_corr_up = 0;
    TH1D *delta_phi_deta3_gap_norm_corr_down = 0;
    TH1D *delta_phi_deta3_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_corr_up);
    if (delta_phi_deta3_gap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta3_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_corr_down);
    if (delta_phi_deta3_gap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta3_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_uncorr);
    if (delta_phi_deta3_gap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta3_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap_norm;
    delta_phi_deta3_gap_norm =  new TH1D("total_unc_delta_phi_deta3_gap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_corr_up, delta_phi_deta3_gap_norm_corr_down, delta_phi_deta3_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down,  8, "top_right", output_path_plots, "tot_unc_delta_phi_deta3_gap_norm", detail);


//compute total uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi Deta4 Gap"<<endl; }

    TH1D *delta_phi_deta4_gap_corr_up = 0;
    TH1D *delta_phi_deta4_gap_corr_down = 0;
    TH1D *delta_phi_deta4_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta4_gap",delta_phi_deta4_gap_corr_up);
    if (delta_phi_deta4_gap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta4_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta4_gap",delta_phi_deta4_gap_corr_down);
    if (delta_phi_deta4_gap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta4_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta4_gap",delta_phi_deta4_gap_uncorr);
    if (delta_phi_deta4_gap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta4_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D("total_unc_delta_phi_deta4_gap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta4_gap, delta_phi_deta4_gap_corr_up, delta_phi_deta4_gap_corr_down, delta_phi_deta4_gap_uncorr, total_unc_up, total_unc_down, 9, "top_right", output_path_plots, "tot_unc_delta_phi_deta4_gap", detail);


//compute total uncertainty for delta phi deta4 gap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Gap Norm"<<endl; }

    TH1D *delta_phi_deta4_gap_norm_corr_up = 0;
    TH1D *delta_phi_deta4_gap_norm_corr_down = 0;
    TH1D *delta_phi_deta4_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_corr_up);
    if (delta_phi_deta4_gap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta4_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_corr_down);
    if (delta_phi_deta4_gap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta4_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_uncorr);
    if (delta_phi_deta4_gap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta4_gap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap_norm;
    delta_phi_deta4_gap_norm =  new TH1D("total_unc_delta_phi_deta4_gap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_corr_up, delta_phi_deta4_gap_norm_corr_down, delta_phi_deta4_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 9, "top_right", output_path_plots, "tot_unc_delta_phi_deta4_gap_norm", detail);


//compute total uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi Nogap"<<endl; }

    TH1D *delta_phi_nogap_corr_up = 0;
    TH1D *delta_phi_nogap_corr_down = 0;
    TH1D *delta_phi_nogap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_nogap",delta_phi_nogap_corr_up);
    if (delta_phi_nogap_corr_up == 0) { cout << "corr_unc_up_delta_phi_nogap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_nogap",delta_phi_nogap_corr_down);
    if (delta_phi_nogap_corr_down == 0) { cout << "corr_unc_down_delta_phi_nogap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_nogap",delta_phi_nogap_uncorr);
    if (delta_phi_nogap_uncorr == 0) { cout << "uncorr_unc_delta_phi_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D("total_unc_delta_phi_nogap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_nogap, delta_phi_nogap_corr_up, delta_phi_nogap_corr_down, delta_phi_nogap_uncorr, total_unc_up, total_unc_down, 10, "top_right", output_path_plots, "tot_unc_delta_phi_nogap", detail);


//compute total uncertainty for delta phi nogap norm distribution
    if (detail) { cout<<"Delta phi Nogap Norm"<<endl; }

    TH1D *delta_phi_nogap_norm_corr_up = 0;
    TH1D *delta_phi_nogap_norm_corr_down = 0;
    TH1D *delta_phi_nogap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_nogap_norm",delta_phi_nogap_norm_corr_up);
    if (delta_phi_nogap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_nogap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_nogap_norm",delta_phi_nogap_norm_corr_down);
    if (delta_phi_nogap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_nogap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_nogap_norm",delta_phi_nogap_norm_uncorr);
    if (delta_phi_nogap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_nogap_norm;
    delta_phi_nogap_norm =  new TH1D("total_unc_delta_phi_nogap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_nogap_norm, delta_phi_nogap_norm_corr_up, delta_phi_nogap_norm_corr_down, delta_phi_nogap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 10, "top_right", output_path_plots, "tot_unc_delta_phi_nogap_norm", detail);


//compute total uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_corr_up = 0;
    TH1D *delta_phi_deta1_nogap_corr_down = 0;
    TH1D *delta_phi_deta1_nogap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta1_nogap",delta_phi_deta1_nogap_corr_up);
    if (delta_phi_deta1_nogap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta1_nogap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta1_nogap",delta_phi_deta1_nogap_corr_down);
    if (delta_phi_deta1_nogap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta1_nogap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta1_nogap",delta_phi_deta1_nogap_uncorr);
    if (delta_phi_deta1_nogap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta1_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D("total_unc_delta_phi_deta1_nogap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta1_nogap, delta_phi_deta1_nogap_corr_up, delta_phi_deta1_nogap_corr_down, delta_phi_deta1_nogap_uncorr, total_unc_up, total_unc_down, 11, "top_right", output_path_plots, "tot_unc_delta_phi_deta1_nogap", detail);


//compute total uncertainty for delta phi deta1 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta1_nogap_norm_corr_up = 0;
    TH1D *delta_phi_deta1_nogap_norm_corr_down = 0;
    TH1D *delta_phi_deta1_nogap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_corr_up);
    if (delta_phi_deta1_nogap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_corr_down);
    if (delta_phi_deta1_nogap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_uncorr);
    if (delta_phi_deta1_nogap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap_norm;
    delta_phi_deta1_nogap_norm =  new TH1D("total_unc_delta_phi_deta1_nogap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_corr_up, delta_phi_deta1_nogap_norm_corr_down, delta_phi_deta1_nogap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 11, "top_right", output_path_plots, "tot_unc_delta_phi_deta1_nogap_norm", detail);


//compute total uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_corr_up = 0;
    TH1D *delta_phi_deta2_nogap_corr_down = 0;
    TH1D *delta_phi_deta2_nogap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta2_nogap",delta_phi_deta2_nogap_corr_up);
    if (delta_phi_deta2_nogap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta2_nogap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta2_nogap",delta_phi_deta2_nogap_corr_down);
    if (delta_phi_deta2_nogap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta2_nogap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta2_nogap",delta_phi_deta2_nogap_uncorr);
    if (delta_phi_deta2_nogap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta2_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D("total_unc_delta_phi_deta2_nogap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta2_nogap, delta_phi_deta2_nogap_corr_up, delta_phi_deta2_nogap_corr_down, delta_phi_deta2_nogap_uncorr, total_unc_up, total_unc_down, 12, "top_right", output_path_plots, "tot_unc_delta_phi_deta2_nogap", detail);


//compute total uncertainty for delta phi deta2 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta2_nogap_norm_corr_up = 0;
    TH1D *delta_phi_deta2_nogap_norm_corr_down = 0;
    TH1D *delta_phi_deta2_nogap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_corr_up);
    if (delta_phi_deta2_nogap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_corr_down);
    if (delta_phi_deta2_nogap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_uncorr);
    if (delta_phi_deta2_nogap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap_norm;
    delta_phi_deta2_nogap_norm =  new TH1D("total_unc_delta_phi_deta2_nogap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_corr_up, delta_phi_deta2_nogap_norm_corr_down, delta_phi_deta2_nogap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 12, "top_right", output_path_plots, "tot_unc_delta_phi_deta2_nogap_norm", detail);


//compute total uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_corr_up = 0;
    TH1D *delta_phi_deta3_nogap_corr_down = 0;
    TH1D *delta_phi_deta3_nogap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta3_nogap",delta_phi_deta3_nogap_corr_up);
    if (delta_phi_deta3_nogap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta3_nogap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta3_nogap",delta_phi_deta3_nogap_corr_down);
    if (delta_phi_deta3_nogap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta3_nogap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta3_nogap",delta_phi_deta3_nogap_uncorr);
    if (delta_phi_deta3_nogap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta3_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D("total_unc_delta_phi_deta3_nogap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta3_nogap, delta_phi_deta3_nogap_corr_up, delta_phi_deta3_nogap_corr_down, delta_phi_deta3_nogap_uncorr, total_unc_up, total_unc_down, 13, "top_right", output_path_plots, "tot_unc_delta_phi_deta3_nogap", detail);


//compute total uncertainty for delta phi deta3 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta3_nogap_norm_corr_up = 0;
    TH1D *delta_phi_deta3_nogap_norm_corr_down = 0;
    TH1D *delta_phi_deta3_nogap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_corr_up);
    if (delta_phi_deta3_nogap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_corr_down);
    if (delta_phi_deta3_nogap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_uncorr);
    if (delta_phi_deta3_nogap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap_norm;
    delta_phi_deta3_nogap_norm =  new TH1D("total_unc_delta_phi_deta3_nogap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_corr_up, delta_phi_deta3_nogap_norm_corr_down, delta_phi_deta3_nogap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 13, "top_right", output_path_plots, "tot_unc_delta_phi_deta3_nogap_norm", detail);


//compute total uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_corr_up = 0;
    TH1D *delta_phi_deta4_nogap_corr_down = 0;
    TH1D *delta_phi_deta4_nogap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta4_nogap",delta_phi_deta4_nogap_corr_up);
    if (delta_phi_deta4_nogap_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta4_nogap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta4_nogap",delta_phi_deta4_nogap_corr_down);
    if (delta_phi_deta4_nogap_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta4_nogap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta4_nogap",delta_phi_deta4_nogap_uncorr);
    if (delta_phi_deta4_nogap_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta4_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D("total_unc_delta_phi_deta4_nogap","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta4_nogap, delta_phi_deta4_nogap_corr_up, delta_phi_deta4_nogap_corr_down, delta_phi_deta4_nogap_uncorr, total_unc_up, total_unc_down, 14, "top_right", output_path_plots, "tot_unc_delta_phi_deta4_nogap", detail);


//compute total uncertainty for delta phi deta4 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta4_nogap_norm_corr_up = 0;
    TH1D *delta_phi_deta4_nogap_norm_corr_down = 0;
    TH1D *delta_phi_deta4_nogap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_corr_up);
    if (delta_phi_deta4_nogap_norm_corr_up == 0) { cout << "corr_unc_up_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_corr_down);
    if (delta_phi_deta4_nogap_norm_corr_down == 0) { cout << "corr_unc_down_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_uncorr);
    if (delta_phi_deta4_nogap_norm_uncorr == 0) { cout << "uncorr_unc_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap_norm;
    delta_phi_deta4_nogap_norm =  new TH1D("total_unc_delta_phi_deta4_nogap_norm","Total Relative Uncertainty;#Delta#phi [rad];Total Relative Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty_v2(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_corr_up, delta_phi_deta4_nogap_norm_corr_down, delta_phi_deta4_nogap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 14, "top_right", output_path_plots, "tot_unc_delta_phi_deta4_nogap_norm", detail);


//compute total uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pT Inside Gap"<<endl; }

    TH1D *leading_pt_inside_gap_corr_up = 0;
    TH1D *leading_pt_inside_gap_corr_down = 0;
    TH1D *leading_pt_inside_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_leading_pt_inside_gap",leading_pt_inside_gap_corr_up);
    if (leading_pt_inside_gap_corr_up == 0) { cout << "corr_unc_up_leading_pt_inside_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_leading_pt_inside_gap",leading_pt_inside_gap_corr_down);
    if (leading_pt_inside_gap_corr_down == 0) { cout << "corr_unc_down_leading_pt_inside_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_leading_pt_inside_gap",leading_pt_inside_gap_uncorr);
    if (leading_pt_inside_gap_uncorr == 0) { cout << "uncorr_unc_leading_pt_inside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D("total_unc_leading_pt_inside_gap","Total Relative Uncertainty;p_{T}^{inside} [GeV];Total Relative Uncertainty", in_nbins, in_bins);

    calc_total_uncertainty_v2(leading_pt_inside_gap, leading_pt_inside_gap_corr_up, leading_pt_inside_gap_corr_down, leading_pt_inside_gap_uncorr, total_unc_up, total_unc_down, 15, "top_right", output_path_plots, "tot_unc_leading_pt_inside_gap", detail);


//compute total uncertainty for leading pt inside gap norm distribution
    if (detail) { cout<<"Leading pT Inside Gap Norm"<<endl; }

    TH1D *leading_pt_inside_gap_norm_corr_up = 0;
    TH1D *leading_pt_inside_gap_norm_corr_down = 0;
    TH1D *leading_pt_inside_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_corr_up);
    if (leading_pt_inside_gap_norm_corr_up == 0) { cout << "corr_unc_up_leading_pt_inside_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_corr_down);
    if (leading_pt_inside_gap_norm_corr_down == 0) { cout << "corr_unc_down_leading_pt_inside_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_uncorr);
    if (leading_pt_inside_gap_norm_uncorr == 0) { cout << "uncorr_unc_leading_pt_inside_gap_norm not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap_norm;
    leading_pt_inside_gap_norm =  new TH1D("total_unc_leading_pt_inside_gap_norm","Total Relative Uncertainty;p_{T}^{inside} [GeV];Total Relative Uncertainty", in_nbins, in_bins);

    calc_total_uncertainty_v2(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_corr_up, leading_pt_inside_gap_norm_corr_down, leading_pt_inside_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 15, "top_right", output_path_plots, "tot_unc_leading_pt_inside_gap_norm", detail);


//compute total uncertainty for leading eta* inside gap distribution
    if (detail) { cout<<"Leading Eta* Gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_corr_up = 0;
    TH1D *leading_eta_star_inside_gap_corr_down = 0;
    TH1D *leading_eta_star_inside_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_leading_eta_star_inside_gap",leading_eta_star_inside_gap_corr_up);
    if (leading_eta_star_inside_gap_corr_up == 0) { cout << "corr_unc_up_leading_eta_star_inside_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_leading_eta_star_inside_gap",leading_eta_star_inside_gap_corr_down);
    if (leading_eta_star_inside_gap_corr_down == 0) { cout << "corr_unc_down_leading_eta_star_inside_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_leading_eta_star_inside_gap",leading_eta_star_inside_gap_uncorr);
    if (leading_eta_star_inside_gap_uncorr == 0) { cout << "uncorr_unc_leading_eta_star_inside_gap not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D("total_unc_leading_eta_star_inside_gap","Total Relative Uncertainty;#eta*;Total Relative Uncertainty", etastar_nbins, etastar_bins);

    calc_total_uncertainty_v2(leading_eta_star_inside_gap, leading_eta_star_inside_gap_corr_up, leading_eta_star_inside_gap_corr_down, leading_eta_star_inside_gap_uncorr, total_unc_up, total_unc_down, 16, "top_right", output_path_plots, "tot_unc_leading_eta_star_inside_gap", detail);


//compute total uncertainty for leading eta* inside gap norm distribution
    if (detail) { cout<<"Leading Eta* Gap Norm"<<endl; }

    TH1D *leading_eta_star_inside_gap_norm_corr_up = 0;
    TH1D *leading_eta_star_inside_gap_norm_corr_down = 0;
    TH1D *leading_eta_star_inside_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_corr_up);
    if (leading_eta_star_inside_gap_norm_corr_up == 0) { cout << "corr_unc_up_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_corr_down);
    if (leading_eta_star_inside_gap_norm_corr_down == 0) { cout << "corr_unc_down_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_uncorr);
    if (leading_eta_star_inside_gap_norm_uncorr == 0) { cout << "uncorr_unc_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap_norm;
    leading_eta_star_inside_gap_norm =  new TH1D("total_unc_leading_eta_star_inside_gap_norm","Total Relative Uncertainty;#eta*;Total Relative Uncertainty", etastar_nbins, etastar_bins);

    calc_total_uncertainty_v2(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_corr_up, leading_eta_star_inside_gap_norm_corr_down, leading_eta_star_inside_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 16, "top_right", output_path_plots, "tot_unc_leading_eta_star_inside_gap_norm", detail);


//compute total uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Delta Eta Outside Gap"<<endl; }

    TH1D *delta_eta_outside_gap_corr_up = 0;
    TH1D *delta_eta_outside_gap_corr_down = 0;
    TH1D *delta_eta_outside_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_eta_outside_gap",delta_eta_outside_gap_corr_up);
    if (delta_eta_outside_gap_corr_up == 0) { cout << "corr_unc_up_delta_eta_outside_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_eta_outside_gap",delta_eta_outside_gap_corr_down);
    if (delta_eta_outside_gap_corr_down == 0) { cout << "corr_unc_down_delta_eta_outside_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_eta_outside_gap",delta_eta_outside_gap_uncorr);
    if (delta_eta_outside_gap_uncorr == 0) { cout << "uncorr_unc_delta_eta_outside_gap not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D("total_unc_delta_eta_outside_gap","Total Relative Uncertainty;#Delta#eta^{out};Total Relative Uncertainty", deta_out_nbins, deta_out_bins);

    calc_total_uncertainty_v2(delta_eta_outside_gap, delta_eta_outside_gap_corr_up, delta_eta_outside_gap_corr_down, delta_eta_outside_gap_uncorr, total_unc_up, total_unc_down, 17, "top_right", output_path_plots, "tot_unc_delta_eta_outside_gap", detail);


//compute total uncertainty for delta eta outside gap norm distribution
    if (detail) { cout<<"Delta Eta Outside Gap Norm"<<endl; }

    TH1D *delta_eta_outside_gap_norm_corr_up = 0;
    TH1D *delta_eta_outside_gap_norm_corr_down = 0;
    TH1D *delta_eta_outside_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_corr_up);
    if (delta_eta_outside_gap_norm_corr_up == 0) { cout << "corr_unc_up_delta_eta_outside_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_corr_down);
    if (delta_eta_outside_gap_norm_corr_down == 0) { cout << "corr_unc_down_delta_eta_outside_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_uncorr);
    if (delta_eta_outside_gap_norm_uncorr == 0) { cout << "uncorr_unc_delta_eta_outside_gap_norm not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap_norm;
    delta_eta_outside_gap_norm =  new TH1D("total_unc_delta_eta_outside_gap_norm","Total Relative Uncertainty;#Delta#eta^{out};Total Relative Uncertainty", deta_out_nbins, deta_out_bins);

    calc_total_uncertainty_v2(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_corr_up, delta_eta_outside_gap_norm_corr_down, delta_eta_outside_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 17, "top_right", output_path_plots, "tot_unc_delta_eta_outside_gap_norm", detail);


//compute total uncertainty for leading pt outside gap distribution
    if (detail) { cout<<"Leading pT Outside Gap"<<endl; }

    TH1D *leading_pt_outside_gap_corr_up = 0;
    TH1D *leading_pt_outside_gap_corr_down = 0;
    TH1D *leading_pt_outside_gap_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_leading_pt_outside_gap",leading_pt_outside_gap_corr_up);
    if (leading_pt_outside_gap_corr_up == 0) { cout << "corr_unc_up_leading_pt_outside_gap not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_leading_pt_outside_gap",leading_pt_outside_gap_corr_down);
    if (leading_pt_outside_gap_corr_down == 0) { cout << "corr_unc_down_leading_pt_outside_gap not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_leading_pt_outside_gap",leading_pt_outside_gap_uncorr);
    if (leading_pt_outside_gap_uncorr == 0) { cout << "uncorr_unc_leading_pt_outside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D("total_unc_leading_pt_outside_gap","Total Relative Uncertainty;p_{T}^{outside} [GeV];Total Relative Uncertainty", out_nbins, out_bins);

    calc_total_uncertainty_v2(leading_pt_outside_gap, leading_pt_outside_gap_corr_up, leading_pt_outside_gap_corr_down, leading_pt_outside_gap_uncorr, total_unc_up, total_unc_down, 18, "top_right", output_path_plots, "tot_unc_leading_pt_outside_gap", detail);


//compute total uncertainty for leading pt outside gap norm distribution
    if (detail) { cout<<"Leading pT Outside Gap Norm"<<endl; }

    TH1D *leading_pt_outside_gap_norm_corr_up = 0;
    TH1D *leading_pt_outside_gap_norm_corr_down = 0;
    TH1D *leading_pt_outside_gap_norm_uncorr = 0;

    data_corr_up->GetObject("corr_unc_up_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_corr_up);
    if (leading_pt_outside_gap_norm_corr_up == 0) { cout << "corr_unc_up_leading_pt_outside_gap_norm not found!" << endl; return; }
    data_corr_down->GetObject("corr_unc_down_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_corr_down);
    if (leading_pt_outside_gap_norm_corr_down == 0) { cout << "corr_unc_down_leading_pt_outside_gap_norm not found!" << endl; return; }
    data_uncorr->GetObject("uncorr_unc_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_uncorr);
    if (leading_pt_outside_gap_norm_uncorr == 0) { cout << "uncorr_unc_leading_pt_outside_gap_norm not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap_norm;
    leading_pt_outside_gap_norm =  new TH1D("total_unc_leading_pt_outside_gap_norm","Total Relative Uncertainty;p_{T}^{outside} [GeV];Total Relative Uncertainty", out_nbins, out_bins);

    calc_total_uncertainty_v2(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_corr_up, leading_pt_outside_gap_norm_corr_down, leading_pt_outside_gap_norm_uncorr, total_unc_norm_up, total_unc_norm_down, 18, "top_right", output_path_plots, "tot_unc_leading_pt_outside_gap_norm", detail);


show_total_uncertainties(total_unc_up);
show_total_uncertainties(total_unc_down);

show_total_uncertainties(total_unc_norm_up);
show_total_uncertainties(total_unc_norm_down);

//results not saved into root file, but they should - also the values are not presented!

}

void total_uncertainty(string path_jes, string path_model, string path_data, string total_uncertainty, string label_in, string label_out, string output_path_plots, bool detail = false, bool disp_uncertainty = true, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Total Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for JES:   " << path_jes << endl; }
    if (detail) { cout << "Input path for Model: " << path_model << endl; }
    if (detail) { cout << "Input path for Stat:  " << path_data << endl; }
    if (detail) { cout << "Label In:             " << label_in << endl; }
    if (detail) { cout << "Label Out:            " << label_out << endl; }
    if (detail) { cout << "Output path:          " << total_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:    " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:         " << detail << endl; }
    if (detail) { cout << "Display Results:      " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:            " << test << endl; }

//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_jes = new TFile( path_jes.c_str() );
    TFile *data_model = new TFile( path_model.c_str() );
    TFile *data = new TFile( path_data.c_str() );


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
    double total_unc[19*3];

    for (int i=0; i<= 19*3-1;i++)
    {
    total_unc[i] = 0.0;
    }


//compute total uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_jes = 0;
    TH1D *delta_phi_model = 0;
    TH1D *delta_phi_data = 0;
    TString delta_phi_name_in = label_in + "delta_phi";
    TString delta_phi_name_out = label_out + "delta_phi";

    data_jes->GetObject(delta_phi_name_in,delta_phi_jes);
    if (delta_phi_jes == 0) { cout << delta_phi_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi",delta_phi_model);
    if (delta_phi_model == 0) { cout << "merged_model_unc_delta_phi not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi",delta_phi_data);
    if (delta_phi_data == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    
    TH1D *delta_phi;
    delta_phi =  new TH1D(delta_phi_name_out,"Total Uncertainty;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi, delta_phi_jes, delta_phi_model, delta_phi_data, total_unc, 0, detail);
    plot_histogram(delta_phi, output_path_plots, label_out + "delta_phi",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_jes = 0;
    TH1D *delta_phi_deta1_model = 0;
    TH1D *delta_phi_deta1_data = 0;
    TString delta_phi_deta1_name_in = label_in + "delta_phi_deta1";
    TString delta_phi_deta1_name_out = label_out + "delta_phi_deta1";

    data_jes->GetObject(delta_phi_deta1_name_in,delta_phi_deta1_jes);
    if (delta_phi_deta1_jes == 0) { cout << delta_phi_deta1_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1",delta_phi_deta1_model);
    if (delta_phi_deta1_model == 0) { cout << "merged_model_unc_delta_phi_deta1 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D(delta_phi_deta1_name_out,"Total Uncertainty deta1;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta1, delta_phi_deta1_jes, delta_phi_deta1_model, delta_phi_deta1_data, total_unc, 1, detail);
    plot_histogram(delta_phi_deta1, output_path_plots, label_out + "delta_phi_deta1",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_jes = 0;
    TH1D *delta_phi_deta2_model = 0;
    TH1D *delta_phi_deta2_data = 0;
    TString delta_phi_deta2_name_in = label_in + "delta_phi_deta2";
    TString delta_phi_deta2_name_out = label_out + "delta_phi_deta2";

    data_jes->GetObject(delta_phi_deta2_name_in,delta_phi_deta2_jes);
    if (delta_phi_deta2_jes == 0) { cout << delta_phi_deta2_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2",delta_phi_deta2_model);
    if (delta_phi_deta2_model == 0) { cout << "merged_model_unc_delta_phi_deta2 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D(delta_phi_deta2_name_out,"Total Uncertainty deta2;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta2, delta_phi_deta2_jes, delta_phi_deta2_model, delta_phi_deta2_data, total_unc, 2, detail);
    plot_histogram(delta_phi_deta2, output_path_plots, label_out + "delta_phi_deta2",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_jes = 0;
    TH1D *delta_phi_deta3_model = 0;
    TH1D *delta_phi_deta3_data = 0;
    TString delta_phi_deta3_name_in = label_in + "delta_phi_deta3";
    TString delta_phi_deta3_name_out = label_out + "delta_phi_deta3";

    data_jes->GetObject(delta_phi_deta3_name_in,delta_phi_deta3_jes);
    if (delta_phi_deta3_jes == 0) { cout << delta_phi_deta3_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3",delta_phi_deta3_model);
    if (delta_phi_deta3_model == 0) { cout << "merged_model_unc_delta_phi_deta3 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D(delta_phi_deta3_name_out,"Total Uncertainty deta3;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta3, delta_phi_deta3_jes, delta_phi_deta3_model, delta_phi_deta3_data, total_unc, 3, detail);
    plot_histogram(delta_phi_deta3, output_path_plots, label_out + "delta_phi_deta3",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_jes = 0;
    TH1D *delta_phi_deta4_model = 0;
    TH1D *delta_phi_deta4_data = 0;
    TString delta_phi_deta4_name_in = label_in + "delta_phi_deta4";
    TString delta_phi_deta4_name_out = label_out + "delta_phi_deta4";

    data_jes->GetObject(delta_phi_deta4_name_in,delta_phi_deta4_jes);
    if (delta_phi_deta4_jes == 0) { cout << delta_phi_deta4_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4",delta_phi_deta4_model);
    if (delta_phi_deta4_model == 0) { cout << "merged_model_unc_delta_phi_deta4 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D(delta_phi_deta4_name_out,"Total Uncertainty deta4;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta4, delta_phi_deta4_jes, delta_phi_deta4_model, delta_phi_deta4_data, total_unc, 4, detail);
    plot_histogram(delta_phi_deta4, output_path_plots, label_out + "delta_phi_deta4",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_jes = 0;
    TH1D *delta_phi_gap_model = 0;
    TH1D *delta_phi_gap_data = 0;
    TString delta_phi_gap_name_in = label_in + "delta_phi_gap";
    TString delta_phi_gap_name_out = label_out + "delta_phi_gap";

    data_jes->GetObject(delta_phi_gap_name_in,delta_phi_gap_jes);
    if (delta_phi_gap_jes == 0) { cout << delta_phi_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_gap",delta_phi_gap_model);
    if (delta_phi_gap_model == 0) { cout << "merged_model_unc_delta_phi_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D(delta_phi_gap_name_out,"Total Uncertainty gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_gap, delta_phi_gap_jes, delta_phi_gap_model, delta_phi_gap_data, total_unc, 5, detail);
    plot_histogram(delta_phi_gap, output_path_plots, label_out + "delta_phi_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_jes = 0;
    TH1D *delta_phi_deta1_gap_model = 0;
    TH1D *delta_phi_deta1_gap_data = 0;
    TString delta_phi_deta1_gap_name_in = label_in + "delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_name_out = label_out + "delta_phi_deta1_gap";

    data_jes->GetObject(delta_phi_deta1_gap_name_in,delta_phi_deta1_gap_jes);
    if (delta_phi_deta1_gap_jes == 0) { cout << delta_phi_deta1_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_gap",delta_phi_deta1_gap_model);
    if (delta_phi_deta1_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta1_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D(delta_phi_deta1_gap_name_out,"Total Uncertainty deta1 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta1_gap, delta_phi_deta1_gap_jes, delta_phi_deta1_gap_model, delta_phi_deta1_gap_data, total_unc, 6, detail);
    plot_histogram(delta_phi_deta1_gap, output_path_plots, label_out + "delta_phi_deta1_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_jes = 0;
    TH1D *delta_phi_deta2_gap_model = 0;
    TH1D *delta_phi_deta2_gap_data = 0;
    TString delta_phi_deta2_gap_name_in = label_in + "delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_name_out = label_out + "delta_phi_deta2_gap";

    data_jes->GetObject(delta_phi_deta2_gap_name_in,delta_phi_deta2_gap_jes);
    if (delta_phi_deta2_gap_jes == 0) { cout << delta_phi_deta2_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_gap",delta_phi_deta2_gap_model);
    if (delta_phi_deta2_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta2_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D(delta_phi_deta2_gap_name_out,"Total Uncertainty deta2 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta2_gap, delta_phi_deta2_gap_jes, delta_phi_deta2_gap_model, delta_phi_deta2_gap_data, total_unc, 7, detail);
    plot_histogram(delta_phi_deta2_gap, output_path_plots, label_out + "delta_phi_deta2_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_jes = 0;
    TH1D *delta_phi_deta3_gap_model = 0;
    TH1D *delta_phi_deta3_gap_data = 0;
    TString delta_phi_deta3_gap_name_in = label_in + "delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_name_out = label_out + "delta_phi_deta3_gap";

    data_jes->GetObject(delta_phi_deta3_gap_name_in,delta_phi_deta3_gap_jes);
    if (delta_phi_deta3_gap_jes == 0) { cout << delta_phi_deta3_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_gap",delta_phi_deta3_gap_model);
    if (delta_phi_deta3_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta3_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D(delta_phi_deta3_gap_name_out,"Total Uncertainty deta3 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta3_gap, delta_phi_deta3_gap_jes, delta_phi_deta3_gap_model, delta_phi_deta3_gap_data, total_unc, 8, detail);
    plot_histogram(delta_phi_deta3_gap, output_path_plots, label_out + "delta_phi_deta3_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_jes = 0;
    TH1D *delta_phi_deta4_gap_model = 0;
    TH1D *delta_phi_deta4_gap_data = 0;
    TString delta_phi_deta4_gap_name_in = label_in + "delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_name_out = label_out + "delta_phi_deta4_gap";

    data_jes->GetObject(delta_phi_deta4_gap_name_in,delta_phi_deta4_gap_jes);
    if (delta_phi_deta4_gap_jes == 0) { cout << delta_phi_deta4_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_gap",delta_phi_deta4_gap_model);
    if (delta_phi_deta4_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta4_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D(delta_phi_deta4_gap_name_out,"Total Uncertainty deta4 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta4_gap, delta_phi_deta4_gap_jes, delta_phi_deta4_gap_model, delta_phi_deta4_gap_data, total_unc, 9, detail);
    plot_histogram(delta_phi_deta4_gap, output_path_plots, label_out + "delta_phi_deta4_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }

    TH1D *delta_phi_nogap_jes = 0;
    TH1D *delta_phi_nogap_model = 0;
    TH1D *delta_phi_nogap_data = 0;
    TString delta_phi_nogap_name_in = label_in + "delta_phi_nogap";
    TString delta_phi_nogap_name_out = label_out + "delta_phi_nogap";

    data_jes->GetObject(delta_phi_nogap_name_in,delta_phi_nogap_jes);
    if (delta_phi_nogap_jes == 0) { cout << delta_phi_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_nogap",delta_phi_nogap_model);
    if (delta_phi_nogap_model == 0) { cout << "merged_model_unc_delta_phi_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D(delta_phi_nogap_name_out,"Total Uncertainty nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_nogap, delta_phi_nogap_jes, delta_phi_nogap_model, delta_phi_nogap_data, total_unc, 10, detail);
    plot_histogram(delta_phi_nogap, output_path_plots, label_out + "delta_phi_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_jes = 0;
    TH1D *delta_phi_deta1_nogap_model = 0;
    TH1D *delta_phi_deta1_nogap_data = 0;
    TString delta_phi_deta1_nogap_name_in = label_in + "delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_name_out = label_out + "delta_phi_deta1_nogap";

    data_jes->GetObject(delta_phi_deta1_nogap_name_in,delta_phi_deta1_nogap_jes);
    if (delta_phi_deta1_nogap_jes == 0) { cout << delta_phi_deta1_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_nogap",delta_phi_deta1_nogap_model);
    if (delta_phi_deta1_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta1_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D(delta_phi_deta1_nogap_name_out,"Total Uncertainty deta1 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jes, delta_phi_deta1_nogap_model, delta_phi_deta1_nogap_data, total_unc, 11, detail);
    plot_histogram(delta_phi_deta1_nogap, output_path_plots, label_out + "delta_phi_deta1_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_jes = 0;
    TH1D *delta_phi_deta2_nogap_model = 0;
    TH1D *delta_phi_deta2_nogap_data = 0;
    TString delta_phi_deta2_nogap_name_in = label_in + "delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_name_out = label_out + "delta_phi_deta2_nogap";

    data_jes->GetObject(delta_phi_deta2_nogap_name_in,delta_phi_deta2_nogap_jes);
    if (delta_phi_deta2_nogap_jes == 0) { cout << delta_phi_deta2_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_nogap",delta_phi_deta2_nogap_model);
    if (delta_phi_deta2_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta2_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D(delta_phi_deta2_nogap_name_out,"Total Uncertainty deta2 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jes, delta_phi_deta2_nogap_model, delta_phi_deta2_nogap_data, total_unc, 12, detail);
    plot_histogram(delta_phi_deta2_nogap, output_path_plots, label_out + "delta_phi_deta2_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_jes = 0;
    TH1D *delta_phi_deta3_nogap_model = 0;
    TH1D *delta_phi_deta3_nogap_data = 0;
    TString delta_phi_deta3_nogap_name_in = label_in + "delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_name_out = label_out + "delta_phi_deta3_nogap";

    data_jes->GetObject(delta_phi_deta3_nogap_name_in,delta_phi_deta3_nogap_jes);
    if (delta_phi_deta3_nogap_jes == 0) { cout << delta_phi_deta3_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_nogap",delta_phi_deta3_nogap_model);
    if (delta_phi_deta3_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta3_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D(delta_phi_deta3_nogap_name_out,"Total Uncertainty deta3 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jes, delta_phi_deta3_nogap_model, delta_phi_deta3_nogap_data, total_unc, 13, detail);
    plot_histogram(delta_phi_deta3_nogap, output_path_plots, label_out + "delta_phi_deta3_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_jes = 0;
    TH1D *delta_phi_deta4_nogap_model = 0;
    TH1D *delta_phi_deta4_nogap_data = 0;
    TString delta_phi_deta4_nogap_name_in = label_in + "delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_name_out = label_out + "delta_phi_deta4_nogap";

    data_jes->GetObject(delta_phi_deta4_nogap_name_in,delta_phi_deta4_nogap_jes);
    if (delta_phi_deta4_nogap_jes == 0) { cout << delta_phi_deta4_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_nogap",delta_phi_deta4_nogap_model);
    if (delta_phi_deta4_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta4_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D(delta_phi_deta4_nogap_name_out,"Total Uncertainty deta4 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jes, delta_phi_deta4_nogap_model, delta_phi_deta4_nogap_data, total_unc, 14, detail);
    plot_histogram(delta_phi_deta4_nogap, output_path_plots, label_out + "delta_phi_deta4_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pt inside gap"<<endl; }

    TH1D *leading_pt_inside_gap_jes = 0;
    TH1D *leading_pt_inside_gap_model = 0;
    TH1D *leading_pt_inside_gap_data = 0;
    TString leading_pt_inside_gap_name_in = label_in + "leading_pt_inside_gap";
    TString leading_pt_inside_gap_name_out = label_out + "leading_pt_inside_gap";

    data_jes->GetObject(leading_pt_inside_gap_name_in,leading_pt_inside_gap_jes);
    if (leading_pt_inside_gap_jes == 0) { cout << leading_pt_inside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_inside_gap",leading_pt_inside_gap_model);
    if (leading_pt_inside_gap_model == 0) { cout << "merged_model_unc_leading_pt_inside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D(leading_pt_inside_gap_name_out,"Total Uncertainty pt leading inside gap;p_{T};Total Uncertainty", in_nbins, in_bins);

    calc_total_uncertainty(leading_pt_inside_gap, leading_pt_inside_gap_jes, leading_pt_inside_gap_model, leading_pt_inside_gap_data, total_unc, 15, detail);
    plot_histogram(leading_pt_inside_gap, output_path_plots, label_out + "leading_pt_inside_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for leading eta star inside gap distribution
    if (detail) { cout<<"Leading eta star inside gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_jes = 0;
    TH1D *leading_eta_star_inside_gap_model = 0;
    TH1D *leading_eta_star_inside_gap_data = 0;
    TString leading_eta_star_inside_gap_name_in = label_in + "leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_name_out = label_out + "leading_eta_star_inside_gap";

    data_jes->GetObject(leading_eta_star_inside_gap_name_in,leading_eta_star_inside_gap_jes);
    if (leading_eta_star_inside_gap_jes == 0) { cout << leading_eta_star_inside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_eta_star_inside_gap",leading_eta_star_inside_gap_model);
    if (leading_eta_star_inside_gap_model == 0) { cout << "merged_model_unc_leading_eta_star_inside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D(leading_eta_star_inside_gap_name_out,"Total Uncertainty eta* leading inside gap;#eta*;Total Uncertainty", etastar_nbins, etastar_bins);

    calc_total_uncertainty(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jes, leading_eta_star_inside_gap_model, leading_eta_star_inside_gap_data, total_unc, 16, detail);
    plot_histogram(leading_eta_star_inside_gap, output_path_plots, label_out + "leading_eta_star_inside_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Delta eta outside gap"<<endl; }

    TH1D *delta_eta_outside_gap_jes = 0;
    TH1D *delta_eta_outside_gap_model = 0;
    TH1D *delta_eta_outside_gap_data = 0;
    TString delta_eta_outside_gap_name_in = label_in + "delta_eta_outside_gap";
    TString delta_eta_outside_gap_name_out = label_out + "delta_eta_outside_gap";

    data_jes->GetObject(delta_eta_outside_gap_name_in,delta_eta_outside_gap_jes);
    if (delta_eta_outside_gap_jes == 0) { cout << delta_eta_outside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_eta_outside_gap",delta_eta_outside_gap_model);
    if (delta_eta_outside_gap_model == 0) { cout << "merged_model_unc_delta_eta_outside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D(delta_eta_outside_gap_name_out,"Total Uncertainty #Delta#eta leading outside gap;#Delta#eta;Total Uncertainty", deta_out_nbins, deta_out_bins);

    calc_total_uncertainty(delta_eta_outside_gap, delta_eta_outside_gap_jes, delta_eta_outside_gap_model, delta_eta_outside_gap_data, total_unc, 17, detail);
    plot_histogram(delta_eta_outside_gap, output_path_plots, label_out + "delta_eta_outside_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for leading pt outside gap distribution
    if (detail) { cout<<"Leading pt outside gap"<<endl; }

    TH1D *leading_pt_outside_gap_jes = 0;
    TH1D *leading_pt_outside_gap_model = 0;
    TH1D *leading_pt_outside_gap_data = 0;
    TString leading_pt_outside_gap_name_in = label_in + "leading_pt_outside_gap";
    TString leading_pt_outside_gap_name_out = label_out + "leading_pt_outside_gap";

    data_jes->GetObject(leading_pt_outside_gap_name_in,leading_pt_outside_gap_jes);
    if (leading_pt_outside_gap_jes == 0) { cout << leading_pt_outside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_outside_gap",leading_pt_outside_gap_model);
    if (leading_pt_outside_gap_model == 0) { cout << "merged_model_unc_leading_pt_outside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D(leading_pt_outside_gap_name_out,"Total Uncertainty pt leading outside gap;p_{T};Total Uncertainty", out_nbins, out_bins);

    calc_total_uncertainty(leading_pt_outside_gap, leading_pt_outside_gap_jes, leading_pt_outside_gap_model, leading_pt_outside_gap_data, total_unc, 18, detail);
    plot_histogram(leading_pt_outside_gap, output_path_plots, label_out + "leading_pt_outside_gap",  label_out + "Uncertainty ", "top_right", true);


//output the error variation
    if (detail) { cout<<"Display the total uncertainties..."<<endl; }
    if (disp_uncertainty) { show_total_uncertainties(total_unc); }


//Opening the output root file
    if (detail) { cout<<"Creating " << total_uncertainty << "..."<<endl; }
    TFile *data_output = TFile::Open( total_uncertainty.c_str() , "RECREATE");

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
    if (detail) { cout<<"Writing was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data_jes->Close();
    data_model->Close();
    data->Close();
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
